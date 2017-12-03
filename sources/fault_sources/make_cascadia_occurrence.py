from make_gsc_oq_inputs import make_collapse_occurrence_text


model = {'src_beta':[-5.0, -6.0, -5.0],
         'beta_wts':[0.68, 0.16, 0.16],
         'src_N0':[0.00230, 0.00230, 0.00230],
         'max_mag':[9.01, 8.98, 9.17],
         'min_mag':8.5,
         'src_weight':1}
         	
mx_dict = {'mx_wts':[0.6, 0.1, 0.3]}
meta = {'one_mx':False, 'beta_wts':[0.68, 0.16, 0.16]}
binwid = 0.1
         	
occtxt = make_collapse_occurrence_text(model, binwid, meta, mx_dict)
print occtxt