# MeasEval → ParamGUI → VICE

ParamGUI R3310 renames the following fields. The preferred source is
`Inv_Param_GUI`, falling back to `Inv_Param` from the complete FluxEval file.
Destinations are under `Parameter`. The GUI itself stays local and is ignored by Git.

| Source | Destination |
| --- | --- |
| `CtrlMaps.Psi_d` | `MACHINE.Ctrl_Map.Psi_d_FIT_param_mat_Vs` |
| `CtrlMaps.Psi_q` | `MACHINE.Ctrl_Map.Psi_q_FIT_param_mat_Vs` |
| `CtrlMaps.L_dd` | `MACHINE.Ctrl_Map.L_dd_param_mat_H` |
| `CtrlMaps.L_qq` | `MACHINE.Ctrl_Map.L_qq_param_mat_H` |
| `CtrlMaps.L_dq` | `MACHINE.Ctrl_Map.L_dq_param_mat_H` |
| `CtrlMaps.L_qd` | `MACHINE.Ctrl_Map.L_qd_param_mat_H` |
| `CtrlMaps.Id_breakpoints` | `MACHINE.Ctrl_Map.Id_fine_param_vec_A` |
| `CtrlMaps.Iq_breakpoints` | `MACHINE.Ctrl_Map.Iq_fine_param_vec_A` |
| `CtrlMaps.EESM_winding_ratio` | `MACHINE.Ctrl_Map.EESM_winding_ratio` |
| `TrqMaps.rpm_breakpoints` | `MACHINE.Trq_Map.rpm_breakpoints_param_vec_rpm` |
| `TrqMaps.Trq_breakpoints` | `MACHINE.Trq_Map.Trq_breakpoints_param_vec_Nm` |
| `TrqMaps.U_DC_breakpoints` | `MACHINE.Trq_Map.U_DC_breakpoints_param_vec_V` |
| `TrqMaps.Id_opti` | `MACHINE.Trq_Map.Id_Opti_param_mat_A` |
| `TrqMaps.Iq_opti` | `MACHINE.Trq_Map.Iq_Opti_param_mat_A` |
| `TrqMaps.I_exc_opti` | `MACHINE.Trq_Map.I_Rot_Opti_param_mat_A` |
| `Limits.Trq_Lim_max` | `MACHINE.Limit_Trq.torque_lim_max_mat_Nm` |
| `Limits.Trq_Lim_min` | `MACHINE.Limit_Trq.torque_lim_min_mat_Nm` |
| `PW.Id_DCI`, `Iq_DCI`, `Id_HFFI`, `Iq_HFFI`, `HFFI_ampl_d`, `HFFI_ampl_q` | Fields with the same names under `PW` |


Torque limits are additionally stored under `MACHINE.Torque_limitation_*`.
In the existing MeasEval export, `Trq_breakpoints` is normalized to −1 through 1,
despite the historical destination suffix `_Nm`.

Maps4D is copied separately from `Inv_Param.Maps4D`:

| Source | Destination |
| --- | --- |
| `CtrlMaps` | `MACHINE.Maps4D.Ctrl_Map` |
| `Interpolated.TrqMaps` | `MACHINE.Maps4D.Trq_Map` |
| `Interpolated.Limits.Trq_Lim_max/min` | `MACHINE.Maps4D.Limit_Trq.torque_lim_max/min_mat_Nm` |
| `Interpolated.PW` | `PW.Maps4D` |


For `MACHINE_GUI.Maps4D`, the save helper uses the GUI variant.

Relevant to the winding factor:

- Legacy EESM path: `Id* = Id + k·Iexc`, `k = EESM_winding_ratio`.
- MATLAB maps: rows are `Iq`, columns are `Id*`. ParamGUI imports copy them.
  Only the VICE Python script transposes the six flux/inductance tables.
- Maps4D uses physical d/q axes for each excitation-current slice. The explicit
  Iexc axis is missing from the inspected `Inv_Param.Maps4D.CtrlMaps` export;
  the original slice order is stored in the LUT.
- GUI saving fills NaNs or replaces them with zeros and masks. It also derives
  `L_dq` from `L_dd` and `L_qd` from `L_qq`.
- The inspected import has no direct assignment from `CtrlMaps.L_exc` to
  `MACHINE.Ctrl_Map.L_Rot_param_mat_H`.
- The separate VICE startup helper mirrors the Id axis. This is a different
  path from the ParamGUI import.

Sources: local ParamGUI R3310, `Functions/MAPS4D_Save.m`, and the Python
parameterization script from VICE `003_004_000_002`. No publication permission
for the complete GUI package was found in the available material, so it remains local.
