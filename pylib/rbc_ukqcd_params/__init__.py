dict_params = {}

import rbc_ukqcd_params.p_test

import rbc_ukqcd_params.p_24D
import rbc_ukqcd_params.p_32D

import rbc_ukqcd_params.p_24DH

import rbc_ukqcd_params.p_16IH2
import rbc_ukqcd_params.p_24IH2

import rbc_ukqcd_params.p_48I

for key, value in dict_params.items():
    value["job_tag"] = key
