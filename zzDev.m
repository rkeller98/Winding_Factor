addpath("./optimization")

eval_set = load_data();
meas = eval_set.measurement;

%%

id = meas.id_req(:);
iq = meas.iq_req(:);
ie = meas.ie_req(:);
psid = meas.psi_d(:);
psiq = meas.psi_q(:);

valid = isfinite(id) & isfinite(iq) & isfinite(ie) & ...
        isfinite(psid) & isfinite(psiq);

id = id(valid);
iq = iq(valid);
ie = ie(valid);
psid = psid(valid);
psiq = psiq(valid);

ie_bp = unique(ie); 

mask_ref = ie == min(ie_bp); 
id_ref = id(mask_ref);
iq_ref = iq(mask_ref);
psi_d_ref = psid(mask_ref);
psi_q_ref = psiq(mask_ref);

id_search = id(~mask_ref);
iq_search = iq(~mask_ref);
ie_search = ie(~mask_ref);
psid_search = psid(~mask_ref);
psiq_search = psiq(~mask_ref);

kd = nan()


for idx = 1:numel(id)
    id_eval = id(idx);
    iq_eval = iq(idx);
    ie_eval = ie(idx);
    psi_d_eval = psid(idx);
    psi_q_eval = psiq(idx);

    
end 
