"""Read grouped Total/Meas data or prepare raw measurements from JSON copies.

Default: Set1's grouped research_flux_measurements.json.gz.
Both importers return the named 1-D arrays in CHANNELS: currents in A,
mechanical speeds in rpm, flux linkages in Vs. Meas fluxes are read directly;
raw samples are averaged and their stationary flux linkages calculated.
run() removes nonfinite measured currents/fluxes with one mask for every array.
Missing optional quantities are NaN and reported together once per import.
save() writes the named arrays to results/python/reduction_data.npz.
"""

import gzip
import json
from pathlib import Path
import warnings

import numpy as np


# Meas channels; reference values belong to the enclosing block.
CHANNELS = {
    "id_meas": "PE_Id_meas_avg",
    "iq_meas": "PE_Iq_meas_avg",
    "ie_meas": "I_exc_meas",
    "rpm_meas": "CAN_RPM_meas",
    "id_req": "PE_Id_req_meas",
    "iq_req": "PE_Iq_req_meas",
    "ie_req": "I_exc_req",
    "rpm_req": "PE_Rpm_req",
    "psi_d": "PE_psi_d_meas",
    "psi_q": "PE_psi_q_meas",
    "ie_ref": "Ref_exc_current",
    "rpm_ref": "Ref_Speed",
}
REQUIRED = ("id_meas", "iq_meas", "ie_meas", "psi_d", "psi_q")


class ReductionDataPreparation:
    def __init__(self, filename=None):
        self.root = Path(__file__).resolve().parents[2]
        self.filename = (
            Path(filename)
            if filename
            else (
                self.root
                / "Resources/DataSets/Set1/extracted/evaluation/research_flux_measurements.json.gz"
            )
        )

    def run(self):
        """Return the twelve named float arrays, filtered together on REQUIRED.

        Missing optional values never remove otherwise valid measurement points.
        """
        with gzip.open(self.filename, "rt") as file:
            loaded = json.load(file)
        if "Total" in loaded:
            data = self.from_meas(loaded["Total"])
        else:
            data = self.from_raw(next(iter(loaded.values())))
        keep = np.isfinite(np.column_stack([data[name] for name in REQUIRED])).all(
            axis=1
        )
        return {name: values[keep] for name, values in data.items()}

    def from_meas(self, total):
        """Read named arrays using CHANNELS; only scalar block references repeat.

        Within each speed group, sort complete blocks stably by Ref_exc_current
        (Resources/calculate_winding_factor.m). Equal references remain separate;
        unknown references follow known ones in input order. Speed-group order
        and measurement order inside each block remain unchanged. rpm_req comes
        from PE_Rpm_req; CAN_RPM_meas is mechanical rpm without conversion.
        Only REQUIRED channels are mandatory. Missing optional quantities become
        NaN, without substituting requested, measured or reference values.
        Length errors identify the block and channel. Filtering happens in run().
        """
        rows = []
        missing = set()
        for speed_name, speed in total.items():
            speed_rows = []
            for block_name, block in speed.items():
                label = f"Total.{speed_name}.{block_name}"
                if "Meas" not in block:
                    raise ValueError(f"{label}: missing Meas channels")
                m = block["Meas"]
                points = {}
                size = None
                for name, channel in CHANNELS.items():
                    reference = name in ("ie_ref", "rpm_ref")
                    source = block if reference else m
                    if channel not in source:
                        if name in REQUIRED:
                            raise ValueError(
                                f"{label}: missing required channel Meas.{channel}"
                            )
                        missing.add(f"{name} ({channel})")
                        column = np.full(size, np.nan)
                    elif reference:
                        value = np.asarray(source[channel], dtype=float)
                        if value.ndim != 0:
                            raise ValueError(
                                f"{label}: {channel} must be a scalar grouping value, got shape {value.shape}"
                            )
                        column = np.full(size, value.item())
                    else:
                        column = np.atleast_1d(np.asarray(source[channel], dtype=float))
                    if column.ndim != 1 or (size is not None and len(column) != size):
                        raise ValueError(
                            f"{label}: {channel} has shape {column.shape}; "
                            f"expected a vector with {size} points as in Meas.PE_Id_meas_avg"
                        )
                    size = len(column)
                    points[name] = column
                current = points["ie_ref"][0] if size else np.nan
                order = (
                    not np.isfinite(current),
                    current if np.isfinite(current) else 0,
                )
                speed_rows.append((order, points))
            rows.extend(
                points for _, points in sorted(speed_rows, key=lambda row: row[0])
            )
        if missing:
            warnings.warn(
                f"{self.filename}: missing optional quantities: {', '.join(sorted(missing))}; using NaN.",
                stacklevel=2,
            )
        return {
            name: np.concatenate([row[name] for row in rows]) if rows else np.empty(0)
            for name in CHANNELS
        }

    def from_raw(self, raw):
        """Average acquisition samples per PE_Rpm_req/I_exc_req/PE_Index_meas.

        Raw current/speed channels are the same as in CHANNELS. Fluxes use
        psi_d=(Uq-Rs*Iq)/omega_el and psi_q=(Rs*Id-Ud)/omega_el, in Vs.
        Raw group references are the recorded I_exc_req/PE_Rpm_req group keys,
        matching matlab/common/prepare_reduction_data.m, not measured currents.
        Keep the existing index>0 and abs(omega_el)>1 preprocessing. Only measured
        currents and flux inputs must be finite before averaging. All other
        channels share those sample selections and arithmetic means.
        Missing/nonfinite group keys leave affected samples unaveraged, with
        unknown references NaN. Within each speed, sort nominal excitation
        groups stably, keeping operating points in first-occurrence order.
        """
        channels = {ch["Name"].replace(" ", ""): ch["Data"] for ch in raw["Y"]}
        names = {
            name: channel
            for name, channel in CHANNELS.items()
            if name not in ("psi_d", "psi_q", "ie_ref", "rpm_ref")
        }
        names.update(
            {
                "ud": "PE_Ud_meas_avg",
                "uq": "PE_Uq_meas_avg",
                "omega": "PE_w_el_rad_avg",
                "resistance": "Const_Machine_R_s",
                "index": "PE_Index_meas",
            }
        )
        required = (
            "id_meas",
            "iq_meas",
            "ie_meas",
            "ud",
            "uq",
            "omega",
            "resistance",
            "index",
        )
        label = f"Raw ({self.filename})"
        data = {}
        missing = set()
        size = None
        for name, channel in names.items():
            if channel not in channels:
                if name in required:
                    raise ValueError(f"{label}: missing required channel {channel}")
                missing.add(f"{name} ({channel})")
                column = np.full(size, np.nan)
            else:
                column = np.atleast_1d(np.asarray(channels[channel], dtype=float))
            if column.ndim != 1 or (size is not None and len(column) != size):
                raise ValueError(
                    f"{label}: {channel} has shape {column.shape}; expected a vector with {size} points as in PE_Id_meas_avg"
                )
            size = len(column)
            data[name] = column

        keep = (data["index"] > 0) & np.isfinite(
            np.column_stack([data[name] for name in required])
        ).all(axis=1)
        values = np.column_stack(list(data.values()))[keep]
        keys = np.column_stack((data["rpm_req"], data["ie_req"], data["index"]))[keep]
        known = np.isfinite(keys[:, :2]).all(axis=1)
        # Unknown groups must not merge unrelated records with the same index.
        keys = np.column_stack(
            (
                np.where(np.isfinite(keys), keys, np.inf),
                np.where(known, -1, np.arange(len(keys))),
            )
        )
        group_keys, first, group = np.unique(
            keys, axis=0, return_index=True, return_inverse=True
        )
        order = np.argsort(group, kind="stable")
        members = (
            np.split(order, np.flatnonzero(np.diff(group[order])) + 1)
            if len(order)
            else []
        )
        averaged = np.array(
            [values[indices].mean(axis=0) for indices in members]
        ).reshape(-1, len(names))
        speed_order = {rpm: i for i, rpm in enumerate(dict.fromkeys(keys[:, 0]))}
        order = sorted(
            range(len(group_keys)),
            key=lambda i: (speed_order[group_keys[i, 0]], group_keys[i, 1], first[i]),
        )
        averaged, group_keys = averaged[order], group_keys[order]
        data = {name: averaged[:, i] for i, name in enumerate(names)}
        data["rpm_ref"] = np.where(
            np.isfinite(group_keys[:, 0]), group_keys[:, 0], np.nan
        )
        data["ie_ref"] = np.where(
            np.isfinite(group_keys[:, 1]), group_keys[:, 1], np.nan
        )
        data = {
            name: column[np.abs(data["omega"]) > 1] for name, column in data.items()
        }
        data["psi_d"] = (data["uq"] - data["resistance"] * data["iq_meas"]) / data[
            "omega"
        ]
        data["psi_q"] = (data["resistance"] * data["id_meas"] - data["ud"]) / data[
            "omega"
        ]
        for name, channel in (("ie_ref", "I_exc_req"), ("rpm_ref", "PE_Rpm_req")):
            if channel not in channels:
                missing.add(f"{name} (group key {channel})")
        if missing:
            warnings.warn(
                f"{self.filename}: missing optional quantities: {', '.join(sorted(missing))}; using NaN.",
                stacklevel=2,
            )
        return {name: data[name] for name in CHANNELS}

    def save(self, data):
        output = self.root / "results/python/reduction_data.npz"
        output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(output, **data)
        return output


if __name__ == "__main__":
    preparation = ReductionDataPreparation()
    data = preparation.run()
    output = preparation.save(data)
    print(f"{len(data['id_meas'])} points saved to {output}")
