from .rbc_ukqcd_params import set_param, get_param

tag = "a_inv_gev"
set_param("test-4nt16", tag, value=1.0)
set_param("24D", tag, value=1.023)
set_param("32D", tag, value=1.023)
set_param("48D", tag, value=1.023)
set_param("48I", tag, value=1.730)
set_param("64I", tag, value=2.359)
set_param("64I-pq", tag, value=2.359)
set_param("32IfineH", tag, value=3.148)

tag = "m_pi"
set_param("test-4nt16", tag, value=0.2)
set_param("24D", tag, value=0.13975)
set_param("32D", tag, value=0.13975)
set_param("48D", tag, value=0.13975)
set_param("48I", tag, value=0.08049)
set_param("32IfineH", tag, value=0.11790)

tag = "m_kk"
set_param("test-4nt16", tag, value=0.4)
set_param("24D", tag, value=0.504154)
set_param("32D", tag, value=0.504154)
set_param("48D", tag, value=0.504154)
set_param("48I", tag, value=0.28853)
set_param("32IfineH", tag, value=0.17720)

tag = "zz_vv"
set_param("test-4nt16", tag, value=0.7)
set_param("24D", tag, value=0.72672)
set_param("32D", tag, value=0.72672)
set_param("48D", tag, value=0.72672)
set_param("48I", tag, value=0.71076)
set_param("64I", tag, value=0.74293)
set_param("64I-pq", tag, value=0.74293)
set_param("32IfineH", tag, value=0.77700)

tag = "zz_aa"
set_param("test-4nt16", tag, value=0.7)
set_param("24D", tag, value=0.73457)
set_param("32D", tag, value=0.73457)
set_param("48D", tag, value=0.73457)
set_param("48I", tag, value=0.71191)
set_param("64I", tag, value=0.74341)
set_param("64I-pq", tag, value=0.74341)
set_param("32IfineH", tag, value=0.77779)

tag = "m_res"
set_param("test-4nt16", tag, value=0.001)
set_param("24D", tag, value=0.0022824)
set_param("32D", tag, value=0.0022824)
set_param("48D", tag, value=0.0022824)
set_param("48I", tag, value=0.0006102)
set_param("64I", tag, value=0.0003116)
set_param("64I-pq", tag, value=0.0003116)
set_param("32IfineH", tag, value=0.0006296)

tag = "m_l"
set_param("test-4nt8", tag, value=0.01)
set_param("test-4nt16", tag, value=0.01)
set_param("24D", tag, value=0.00107)
set_param("32D", tag, value=0.00107)
set_param("48D", tag, value=0.00107)
set_param("48I", tag, value=0.00078)
set_param("64I", tag, value=0.000678)
set_param("64I-pq", tag, value=0.0006203)
set_param("32IfineH", tag, value=0.0047)

tag = "m_h"
set_param("test-4nt8", tag, value=0.04)
set_param("test-4nt16", tag, value=0.04)
set_param("24D", tag, value=0.0850)
set_param("32D", tag, value=0.0850)
set_param("48D", tag, value=0.0850)
set_param("48I", tag, value=0.0362)
set_param("64I", tag, value=0.02661)
set_param("64I-pq", tag, value=0.02539)
set_param("32IfineH", tag, value=0.0186)

# PHYSICAL REVIEW D 93, 074505 (2016)
# zz_m_l * m_l => m_l in MSbar scheme 3 GeV
tag = "zz_m_l"
set_param("64I", tag, value=2.997 / 2.198)
set_param("64I-pq", tag, value=2.997 / 2.198)
set_param("48I", tag, value=2.997 / 2.198 * 0.9715)

# PHYSICAL REVIEW D 93, 074505 (2016)
# zz_m_h * m_h => m_l in MSbar scheme 3 GeV
tag = "zz_m_h"
set_param("64I", tag, value=81.64 / 60.62)
set_param("64I-pq", tag, value=81.64 / 60.62)
set_param("48I", tag, value=81.64 / 60.62 * 0.9628)

tag = "zz_ss_l"
set_param("64I", tag, value=1 / get_param("64I", "zz_m_l"))
set_param("64I-pq", tag, value=1 / get_param("64I-pq", "zz_m_l"))
set_param("48I", tag, value=1 / get_param("48I", "zz_m_l"))

tag = "zz_ss_h"
set_param("64I", tag, value=1 / get_param("64I", "zz_m_h"))
set_param("64I-pq", tag, value=1 / get_param("64I-pq", "zz_m_h"))
set_param("48I", tag, value=1 / get_param("48I", "zz_m_h"))
