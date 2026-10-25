# Measurement data

| Set | Speed | Excitation current | Evaluation |
| --- | --- | --- | --- |
| Set1 | 3000 rpm | 2–40 A | 2 variants available |
| Set2 | 3000 rpm | 0–40 A | 3 variants available |
| Set3 | 2250 rpm | 3–30 A | New research evaluation available |

- `raw/`: unchanged original files.
- `FluxEval/`: complete evaluations; kept locally and ignored by Git.
- `extracted/raw/`: all raw channels grouped by speed and excitation current.
- `extracted/evaluation/`: all saved `Meas` fields from the evaluations.

All raw files and extracted measurements are also available as gzip-compressed
UTF-8 JSON (`.json.gz`) with the same base filename. MAT files are retained unchanged.
JSON objects correspond to the variables and fields in the MAT files; nothing
is filtered or averaged. Nonfinite values (`NaN`, `Inf`) are stored as `null`.
Vectors become lists and matrices become lists of rows; JSON does not retain a
vector's original row/column orientation. `SampleIndex` remains MATLAB-based
(first sample = 1).

Read in Python without additional packages:

```python
import gzip
import json
from pathlib import Path

path = Path("Resources/DataSets/Set1/extracted/evaluation/research_flux_measurements.json.gz")
with gzip.open(path, "rt", encoding="utf-8") as file:
    data = json.load(file)

group = next(iter(next(iter(data["Total"].values())).values()))
measurements = group["Meas"]
```

For raw files, the original measurement structure is under its original variable
name (`next(iter(data.values()))`); `Y` contains the channels with `Name`, `Unit`,
and `Data`, and `X` contains the time axis. For grouped measurements, `Total`
provides direct access without MATLAB cell indexing of `Meas`.

Regenerate JSON copies in MATLAB:

```matlab
addpath('Resources/DataSets');
export_measurements_json;
```

`matlab/common/prepare_reduction_data.m` rebuilds `research_flux_measurements.mat` for
all three sets from the raw channels. It contains arithmetically averaged
operating points, measured currents, and stationary flux linkages calculated
from measured dq voltages. It is a compact research export without production
GUI or LUT state. Set3 is therefore included in the current study. Other
historical exports are preserved.

Extract a single FluxEval file as MAT:

```matlab
addpath('C:\Git\Winding_Factor\Resources');
measFile = extract_meas_from_fluxeval('C:\path\FluxEval.mat');
data = load(measFile);
id = data.Meas.PE_Id_meas_avg{1,1};
```

Optional second argument: the desired output MAT path. For sources in
`FluxEval/`, output is written to `extracted/evaluation/`. Otherwise,
`<original_name>_measurements.mat` is created alongside the source.

`Meas.<channel>{excitation_row,speed_column}` contains all fields in the cell
layout used by `calculate_winding_factor.m`. Speed and excitation current are
sorted numerically; `rpm_ref` and `exc_ref` contain their reference values.
`Total.<rpm>.<Iexc>.Meas` is also retained, including matrices and outliers.

Extract all local files again:

```matlab
addpath('C:\Git\Winding_Factor\Resources\DataSets');
extract_measurements;
```

Raw data retain all repetitions and channels: 268 in Set1/2 and 195 in Set3.
`SampleIndex` and `Time` refer to the original samples. Nothing is averaged or
filtered. The 284 index-0 samples in Set3 are also retained. Existing evaluations
each contain 303 `Meas` fields.
