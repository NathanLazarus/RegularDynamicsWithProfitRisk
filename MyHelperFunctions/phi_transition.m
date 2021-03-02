function phi_plus = phi_transition(phi, phi_lag, Pp, P, LAMBDAphi, LAMBDAphi_lagged, phi_dependence_on_delta_theta)
theta_plus = (Pp/(Pp - 1));
theta_now = (P/(P - 1));
phi_plus = phi ^ LAMBDAphi * phi_lag ^ LAMBDAphi_lagged * (theta_plus - theta_now) ^ phi_dependence_on_delta_theta;