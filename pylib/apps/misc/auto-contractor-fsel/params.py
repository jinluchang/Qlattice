import rbc_ukqcd_params as rup

tag = "n_exact_wsrc"
rup.dict_params["test-4nt8"][tag] = 2
rup.dict_params["48I"][tag] = 2

tag = "prob_exact_wsrc"
rup.dict_params["test-4nt16"][tag] = 1/8
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["32IfineH"][tag] = 1/32
rup.dict_params["24IH1"][tag] = 1/32
rup.dict_params["24IH2"][tag] = 1/32
rup.dict_params["32IH2"][tag] = 1/32

tag = "n_rand_u1_fsel"
rup.dict_params["test-4nt8"][tag] = 4
rup.dict_params["test-4nt16"][tag] = 4
rup.dict_params["48I"][tag] = 16
rup.dict_params["64I"][tag] = 16
rup.dict_params["16IH2"][tag] = 16
rup.dict_params["32IfineH"][tag] = 64
rup.dict_params["24IH1"][tag] = 64
rup.dict_params["24IH2"][tag] = 64
rup.dict_params["32IH2"][tag] = 64

tag = "prob_acc_1_rand_u1"
rup.dict_params["test-4nt8"][tag] = 1/4
rup.dict_params["test-4nt16"][tag] = 1/4
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["32IfineH"][tag] = 1/32
rup.dict_params["24IH1"][tag] = 1/32
rup.dict_params["24IH2"][tag] = 1/32
rup.dict_params["32IH2"][tag] = 1/32

tag = "prob_acc_2_rand_u1"
rup.dict_params["test-4nt8"][tag] = 1/16
rup.dict_params["test-4nt16"][tag] = 1/16
rup.dict_params["16IH2"][tag] = 1/64
rup.dict_params["32IfineH"][tag] = 1/128
rup.dict_params["24IH1"][tag] = 1/128
rup.dict_params["24IH2"][tag] = 1/128
rup.dict_params["32IH2"][tag] = 1/128

tag = "trajs"
rup.dict_params["test-4nt8"][tag] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"][tag] = list(range(1000, 1400, 100))
rup.dict_params["48I"][tag] = list(range(3000, 500, -5))
rup.dict_params["24D"][tag] = list(range(1000, 10000, 10))
rup.dict_params["16IH2"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IfineH"][tag] = list(range(1000, 10000, 50))
rup.dict_params["24IH2"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24IH1"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IH2"][tag] = list(range(1000, 10000, 10))

