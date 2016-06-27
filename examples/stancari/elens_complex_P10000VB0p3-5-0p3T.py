# AUTHOR: Vince Moens
# PROJECT: Master Thesis EPFL 2013

# NOTES ON EMITTANCE TYPES
# GUN: The code automatically injects particles from the cathode conductor using Child Langmuir law. 
# PROFILE: Profiles measured in the test bench are injected into the lattice. The gun is ommitted fromt the lattice. Injection takes place at the end of that anode. 

############################
# >>>  Package Loading <<< #
############################

from warp import * # warp code
from datetime import *
import numpy as np

#########################
# >>>  File Loading <<< #
#########################

# --- Profiles --- #
fns = [
# new profiles acquired with ACL script
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121218_9p25A_3-3-3kG_500V_51mA_b_57_8102_particles.txt", #0
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121218_9p25A_3-3-3kG_8kV_2940mA_58_8099_particles.txt", #1
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121219_9p25A_3-3-3kG_4kV_1100mA_59_8086_particles.txt", #2
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121219_9p25A_3-3-3kG_2kV_388mA_60_8019_particles.txt", #3
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121219_9p25A_3-3-3kG_1kV_140mA_61_8054_particles.txt", #4
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121219_9p25A_3-3-3kG_6kV_1924mA_62_8087_particles.txt", #5
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121219_9p25A_3-3-3kG_3kV_694mA_63_8035_particles.txt", #6
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121219_9p25A_3-3-3kG_7kV_2394mA_64_8082_particles.txt", #7
# tem/../HG1bpPrature-limited regime
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121220_7p25A_3-3-3kG_8kV_440mA_coarse_65_8183_particles.txt", #8
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_121220_7p25A_3-3-3kG_8kV_27mA_66_7901_particles.txt", #9
# spa/../HG1bcP-charge limited regime
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130220_9p25A_1-1-1kG_889V_118mA_67_8177_particles.txt", #10
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130221_9p25A_1-1-1kG_2kV_398mA_68_8071_particles.txt", #11
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130221_9p25A_1-1-1kG_3kV_726mA_69_8041_particles.txt", #12
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130404_9p25A_25-25-25kG_5556V_1770mA_71_8076_particles.txt", #13
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130221_9p25A_2-2-2kG_3556V_928mA_70_8137_particles.txt", #14
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130404_9p25A_25-25-25kG_5kV_1484mA_72_8175_particles.txt", #15
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130404_9p25A_25-25-25kG_2083V_412mA_73_8059_particles.txt", #16
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130404_9p25A_25-25-25kG_3472V_865mA_74_8049_particles.txt", #17
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130404_9p25A_3-3-3kG_5kV_1460mA_75_8156_particles.txt", #18
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130405_9p25A_2-2-2kG_1333V_213mA_76_8072_particles.txt", #19
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130405_9p25A_2-2-2kG_2222V_448mA_77_8027_particles.txt", #20
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130405_9p25A_2-2-2kG_3111V_744mA_78_8050_particles.txt", #21
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130405_9p25A_15-15-15kG_750V_90mA_79_8081_particles.txt", #22
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130405_9p25A_15-15-15kG_2kV_390mA_80_8089_particles.txt", #23
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130405_9p25A_15-15-15kG_1250V_192mA_81_8045_particles.txt", #24
  "../../HG1b/Profile/Results/Chart_Colors/G_091027_775A_3kV_303030kG_58_06_7927_particles.txt", #25
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130502_9p25A_25-25-25kG_500V_51mA_combi_84_8192_particles.txt", #26
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130507_9p25A_3-3-3kG_500V_5_96_8151_particles.txt", #27
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130513_10A_3-3-3kG_500V_50mA_97_8184_particles.txt", #28
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130513_10A_3-3-3kG_500V_50mA_300Hz_98_8206_particles.txt", #29
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130513_10A_3-3-3kG_500V_50mA_99_8184_particles.txt", #30
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130513_10A_3-3-3kG_500V_50mA_300Hz_100_8206_particles.txt", #31
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130513_9p25A_1-3-3kG_500V_70mA_coarse_101_8192_particles.txt", #32
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130513_9p25A_1-3-2kG_500V_70mA_coarse_102_8190_particles.txt", #33
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130514_9p25A_08-32-08kG_500V_73mA_104_8121_particles.txt", #34
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130514_9p25A_08-32-08kG_2kV_526mA_106_8140_particles.txt", #35
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130514_9p25A_08-32-08kG_4kV_1450mA_107_8157_particles.txt", #36
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130514_9p25A_08-32-08kG_1kV_195mA_105_8136_particles.txt", #37
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_08-32-08kG_8kV_3700mA_108_8125_particles.txt", #38
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_08-32-08kG_7kV_3100mA_109_8121_particles.txt", #39
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_08-32-08kG_250V_26mA_110_8122_particles.txt", #40
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_08-32-08kG_125V_9mA_111_8096_particles.txt", #41
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_08-32-08kG_3kV_944mA_112_8161_particles.txt", #42
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_08-32-08kG_5kV_1990mA_113_8114_particles.txt", #43
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_08-32-08kG_6kV_2530mA_114_8113_particles.txt", #44
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_1-4-1kG_3125V_1022mA_116_8113_particles.txt", #45
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_1-4-1kG_1kV_196mA_115_8074_particles.txt", #46
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130515_9p25A_1-4-1kG_6250V_2754mA_117_8192_particles.txt", #47
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130521_9p25A_06-24-06kG_500V_73mA_133_8117_particles.txt", #48
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130516_9p25A_04-16-04kG_1kV_195mA_122_8177_particles.txt", #49
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130516_9p25A_04-16-04kG_500V_72mA_123_8157_particles.txt", #50
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130520_9p25A_1-4-1kG_500V_73mA_124_500_particles.txt", #51
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130520_9p25A_1-4-1kG_2kV_534mA_125_8136_particles.txt", #52
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130520_9p25A_1-4-1kG_4kV_1440mA_126_8115_particles.txt", #53
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130516_9p25A_06-24-06kG_1kV_196mA_118_8075_particles.txt", #54
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130520_9p25A_1-4-1kG_5kV_1980mA_127_8108_particles.txt", #55
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130520_9p25A_04-16-04kG_2kV_488mA_129_8099_particles.txt", #56
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130520_9p25A_1-4-1kG_8kV_3880mA_128_8106_particles.txt", #57
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130520_9p25A_04-16-04kG_750V_125mA_130_8169_particles.txt", #58
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130521_9p25A_04-16-04kG_1500V_332mA_131_8192_particles.txt", #59
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130521_9p25A_04-16-04kG_2500V_658mA_132_8134_particles.txt", #60
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130521_9p25A_06-24-06kG_2kV_526mA_134_8147_particles.txt", #61
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130521_9p25A_06-24-06kG_3kV_924mA_135_8129_particles.txt", #62
  "../../HG1b/Profile/Results/Chart_Colors/HG1b_130521_9p25A_06-24-06kG_4kV_1368mA_136_8192_particles.txt" #63
]

#Bgun=[3,3,3,3,3,3,3,3,3,3,1,1,1,2.5,2,2.5,2.5,2.5,2.5,2,2,2,1.5,1.5,1.5,3,2.5,3,3,3,3,3,1,1,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,1,1,1,0.6,0.4,0.4,1,1,1,0.6,1,0.4,1,0.4,0.4,0.4,0.6,0.6,0.6]
#Bmain=[3,3,3,3,3,3,3,3,3,3,1,1,1,2.5,2,2.5,2.5,2.5,2.5,2,2,2,1.5,1.5,1.5,3,2.5,3,3,3,3,3,3,3,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,4,4,4,2.4,1.6,1.6,4,4,4,2.4,4,1.6,4,1.6,1.6,1.6,2.4,2.4,2.4]
#Bcoll=[3,3,3,3,3,3,3,3,3,3,1,1,1,2.5,2,2.5,2.5,2.5,2.5,2,2,2,1.5,1.5,1.5,3,2.5,3,3,3,3,3,3,2,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,1,1,1,0.6,0.4,0.4,1,1,1,0.6,1,0.4,1,0.4,0.4,0.4,0.6,0.6,0.6]

#Voltage=[500,8000,4000,2000,1000,6000,3000,7000,8000,8000,889,2000,3000,5556,3556,5000,2083,3472,5000,1333,2222,3111,750,2000,1250
#          ,3000,500,500,500,500,500,500,500,500,500,2000,4000,1000,8000,7000,250,125,3000,5000,6000,3125,1000,6250,500,1000,500,500,2000,4000,1000,5000,2000
#          ,8000,750,1500,2500,2000,3000,4000]
Voltage=[500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000]
#Current=[0.051,2.940,1.100,0.388,0.140,1.924,0.694,2.394,0.440,0.027,0.118,0.398,0.726,1.770,0.928,1.484,0.412,0.865,1.460,0.213
#          ,0.448,0.744,0.090,0.390,0.192,0.582,0.051,0.050,0.050,0.050,0.050,0.050,0.070,0.070,0.073,0.526,1.450,0.195,3.7,3.1,0.026,
#          0.009,0.944,1.99,2.53,1.022,0.196,2.754,0.073,0.195,0.072,0.073,0.534,1.44,0.196,1.98,0.488,3.88,0.125,0.332,0.658,0.526,0.924,1.368]

# --- Selecting profile --- #
item = 10
file_ending = "test_v"
Bmain = 5 #float(Bmain[item])/10
print("\nMagetic Field in Main Solenoid: %g T" % Bmain)
Bgun = 0.3 #float(Bgun[item])/10
Bcoll = 0.3 #float(Bcoll[item])/10
Bbend = Bcoll/9.53
Cathode_Potential = -Voltage[item]  #e-3*kV
print("Cathode Potential: %g V" % Cathode_Potential)
#Current=Current[item]
Current= pow(Voltage[item],1.5)*5.3e-6 #Current[item]
print("Current: %g A (approximate, exact value determined by CLL" % Current)
npart= 500*Current/0.06 
compact_factor = 20
file_ending = "P"+str(-Cathode_Potential)+"VB"+str(int(Bgun*10))+"-"+str(int(Bmain*10))+"-"+str(int(Bcoll*10))+"kG"

# if (round(Bmain,1) == 0.4):
#   Xoff = -3.208143e-3
#   Yoff = -1.03759e-3
# elif (round(Bmain,2) == 0.32):
#   Xoff = -3.072298e-3
#   Yoff = -0.926263e-3
# elif (round(Bmain,2) == 0.24):
#   Xoff = -2.381263e-3
#   Yoff = -0.6774335e-3
# elif (round(Bmain,2) == 0.16):
#   Xoff = -1.210455e-3
#   Yoff = -0.4880404e-3
# else:
#    Xoff = 0.0
#    Yoff = 0.0

# print ("X-offset is "+str(Xoff))
# print ("Y-offset is "+str(Yoff))

###################
# >>> Options <<< #
###################

machine_type            = "TEL2s"     # Options are tbench, TEL2s and TEL2
print("Machine type is: "+machine_type+" setup")
if ((machine_type != "tbench") and (machine_type != "TEL2") and (machine_type != "TEL2s")):
  print("Wrong machine type!")
  quit()
machine_injtype        = "gun"        # Options: "profile" or "gun"
print("Injection type is: "+machine_injtype+" injection \n")
if (machine_injtype=="gun"):
  machine_emittype         = 2          # 2: space-charge limited (Child-Langmuir)
elif (machine_injtype=="profile"):
  machine_emittype         = 1          # 1: constant current,
else: 
  print("Wrong machine_injtype!!")
  quit()

####################
# >>>  Headers <<< #
####################

now=datetime.now()
date=now.strftime("%y%m%d")
time=now.strftime("%H%M")


#if not os.path.exists("../Results/"+date+"/"):
#    os.makedirs("../Results/"+date+"/")
#    print("New day folder created \n")
os.system("cp elens_complex_P"+str(-Cathode_Potential)+"VB"+str(int(Bgun*10))+"-"+str(int(Bmain*10))+"-"+str(int(Bcoll*10))+"kG.py ../Results/"+date+"/"+machine_type+"_"+date+time+"_P"+str(-Cathode_Potential)+"VB"+str(int(Bgun*10))+"-"+str(int(Bmain*10))+"-"+str(int(Bcoll*10))+"kG.py")
os.system("cp elens_complex_P"+str(-Cathode_Potential)+"VB"+str(int(Bgun*10))+"-"+str(int(Bmain*10))+"-"+str(int(Bcoll*10))+"kG.run ../Results/"+date+"/"+machine_type+"_"+date+time+"_P"+str(-Cathode_Potential)+"VB"+str(int(Bgun*10))+"-"+str(int(Bmain*10))+"-"+str(int(Bcoll*10))+"kG.run")


top.runid               = machine_type+"_"+date+time+"_"+machine_injtype+"_"+file_ending
if machine_type == "tbench":
  top.pline2            = "Electron Lens Test Bench"
else: top.pline2        = "Tevatron Electron Lens 2"
if machine_emittype == 1:
  top.pline1            = "Constant-injection_" + machine_injtype
elif machine_emittype == 2: 
  top.pline1            = "Child-Langmuir_" + machine_injtype
else: top.pline1        = "other injection method"
top.runmaker            = "V. Moens"


#####################
# >>> Variables <<< #
#####################

# --- Machine Parameters --- #
machine_zstart          = .0e0
if machine_type   == "tbench":
  machine_syslen        = 2.86                      # from anode to collector [m]
elif ((machine_type == "TEL2") or (machine_type == "TEL2s")):
  machine_syslen        = 4.68581
machine_zplat           = machine_syslen            # diagnostic screen [m]
machine_piperad         = 3*cm                     # inner pipe radius [m]
zfinal                  = machine_zstart + machine_syslen           # Position of collector of pinnhole [m]

# --- Electron Gun --- #
# - Cathode
Cathode_zstart          = -29.25*mm
Cathode_zend            = 0.0*mm
Cathode_radi            = 6.75*mm
Cathode_rado            = 12.7*mm
Cathode_radcurvb        = 10*mm
Cahtode_radcurvs        = 0.5*mm
Cathode_voltage         = Cathode_Potential                        # Cathode Votlage [V]
# - Anode
#The values were taken from a drawing printed on tabloid paper and a conversion rate of 2.25mm(real)/mm(drawing)
Anode_zstart            = 9.48*mm
Anode_z1                = Anode_zstart + 1.5*mm
Anode_z2                = Anode_zstart + 3.5*mm
Anode_z3                = Anode_z1 + 9*mm
Anode_z4                = Anode_z3 + 11.25*mm
Anode_z5                = Anode_z4 + 5.625*mm
Anode_zend              = Anode_z5 + 58.5*mm
Anode_ri                = 14.25*mm
Anode_ro                = Anode_ri+5.33*mm
Anode_r1                = Anode_ri
Anode_r2                = Anode_ro
Anode_r3                = Anode_ri
Anode_r4                = Anode_ri + 0.675*mm
Anode_radtipi           = Anode_ri + 1.5*mm
Anode_radtipo           = Anode_ro - 3.5*mm
Anode_r5                = Anode_ri + 5.625*mm
Anode_rendi             = Anode_r5
Anode_rendo             = Anode_rendi + 1.35*mm
Anode_radcurvb          = 3.5*mm
Anode_radcurvs          = -1.5*mm
Anode_voltage           = 0.0e0
# - Electrode F
ElectrodeF_zstart       = Cathode_zstart
ElectrodeF_zend         = 0.98*mm
ElectrodeF_z1           = ElectrodeF_zend - 0.5*mm
ElectrodeF_z2           = ElectrodeF_zend - 1.4*mm
ElectrodeF_ri           = 13.1*mm
ElectrodeF_ro           = ElectrodeF_ri + 1.9*mm
ElectrodeF_r1           = ElectrodeF_ri + 0.5*mm
ElectrodeF_radcurvs     = -0.5*mm
ElectrodeF_radcurvb     = 1.4*mm
ElectrodeF_voltage      = Cathode_Potential
# - Electrode C
ElectrodeC_zstart       = Cathode_zstart
ElectrodeC_zend         = 1.97*mm
ElectrodeC_ri           = 20.5*mm
ElectrodeC_ro           = 22.0*mm
ElectrodeC_radcurv      = 0.75*mm
ElectrodeC_z1           = ElectrodeC_zend-0.75*mm
ElectrodeC_voltage      = Cathode_Potential
# - Gun drift pipe
Gun_pipe_zstart         = 84.375*mm
Gun_pipe_zend           = 178.875*mm
Gun_pipe_ri             = 36*mm  #Should this not be 3 cm?
Gun_pipe_ro             = 33.75*mm
Gun_pipe_voltage        = 0.0

if machine_type == "tbench":
  # ---- Solenoids --- #
  # - Gun Solenoid
  tbench_solenoid_gun_zstart  = -13*cm
  tbench_solenoid_gun_zend    = 37*cm
  tbench_solenoid_gun_radi    = 28*cm
  tbench_solenoid_gun_rado    = tbench_solenoid_gun_radi+5.433*cm
  tbench_solenoid_gun_b       = Bgun                      # Maximum axial B field [T]

  # - Main Solenoid
  tbench_solenoid_main_zstart = 0.60
  tbench_solenoid_main_zend   = 2.52
  tbench_solenoid_main_radi   = 0.20
  tbench_solenoid_main_rado   = tbench_solenoid_main_radi+14.48*cm
  tbench_solenoid_main_b      = Bmain                      # Maximum axial B field [T]

  # - Collector Solenoid
  tbench_solenoid_col_zstart  = 2.67
  tbench_solenoid_col_zend    = 3.17
  tbench_solenoid_col_radi    = 28*cm
  tbench_solenoid_col_rado    = tbench_solenoid_col_radi + 5.433*cm
  tbench_solenoid_col_b       = Bcoll                       # Maximum axial B field [T]

  # --- Drift Spaces --- #
  # - First Drift
  tbench_drift1_zstart        = 37*cm
  tbench_drift1_zend          = 0.60
  tbench_drift1_ap            = machine_piperad 
  # - Second Drift
  tbench_drift2_zstart        = 2.52
  tbench_drift2_zend          = 2.67
  tbench_drift2_ap            = machine_piperad 

elif ((machine_type == "TEL2") or (machine_type == "TEL2s")):
  # --- Solenoids --- #
  # - Gun Solenoid
  TEL2_solenoid_gun_zstart    = -167.1*mm
  TEL2_solenoid_gun_length    = 330*mm
  TEL2_solenoid_gun_zend      = TEL2_solenoid_gun_zstart+TEL2_solenoid_gun_length
  TEL2_solenoid_gun_radi      = 120*mm
  TEL2_solenoid_gun_rado      = 248*mm
  TEL2_solenoid_gun_b         = Bgun
  
  # --- Bend Solenoids --- #
  # - first bend starting from gun
  TEL2_bendsol1_gun_zstart     = Cathode_zend + 281.6*mm
  TEL2_bendsol1_gun_zend       = TEL2_bendsol1_gun_zstart + 90*mm
  TEL2_bendsol1_gun_length     = 90*mm
  TEL2_bendsol1_gun_ri         = 193*mm
  TEL2_bendsol1_gun_ro         = 265*mm
  TEL2_bendsol1_gun_b          = Bbend

  # - first bend starting from gun
  TEL2_bendsol2_gun_zstart     = TEL2_bendsol1_gun_zend + 52.9*mm
  TEL2_bendsol2_gun_zend       = TEL2_bendsol2_gun_zstart + 90*mm
  TEL2_bendsol2_gun_length     = 90*mm
  TEL2_bendsol2_gun_ri         = 193*mm
  TEL2_bendsol2_gun_ro         = 265*mm
  TEL2_bendsol2_gun_b          = Bbend

  # - first bend starting from gun
  TEL2_bendsol3_gun_zstart     = TEL2_bendsol2_gun_zend + 52.9*mm
  TEL2_bendsol3_gun_zend       = TEL2_bendsol3_gun_zstart + 90*mm
  TEL2_bendsol3_gun_length     = 90*mm
  TEL2_bendsol3_gun_ri         = 193*mm
  TEL2_bendsol3_gun_ro         = 265*mm
  TEL2_bendsol3_gun_b          = Bbend

  # - Main Solenoid
  TEL2_solenoid_main_zstart    = TEL2_bendsol3_gun_zend + 82.8*mm
  TEL2_solenoid_main_length    = 2688.5*mm
  TEL2_solenoid_main_zend      = TEL2_solenoid_main_zstart+TEL2_solenoid_main_length
  TEL2_solenoid_main_radi      = 42.75*mm
  TEL2_solenoid_main_rado      = 241*mm
  TEL2_solenoid_main_b         = Bmain

  # --- Bend Solenoids --- #            ### --- These are in a linear alignment. We need to put them in a bent alignment
  # - first bend starting from gun
  TEL2_bendsol1_col_zstart     = TEL2_solenoid_main_zend + 82.8*mm
  TEL2_bendsol1_col_zend       = TEL2_bendsol1_col_zstart + 90*mm
  TEL2_bendsol1_col_length     = 90*mm
  TEL2_bendsol1_col_ri         = 193*mm
  TEL2_bendsol1_col_ro         = 265*mm
  TEL2_bendsol1_col_b          = Bbend

  # - first bend starting from gun
  TEL2_bendsol2_col_zstart     = TEL2_bendsol1_col_zend + 52.9*mm
  TEL2_bendsol2_col_zend       = TEL2_bendsol2_col_zstart + 90*mm
  TEL2_bendsol2_col_length     = 90*mm
  TEL2_bendsol2_col_ri         = 193*mm
  TEL2_bendsol2_col_ro         = 265*mm
  TEL2_bendsol2_col_b          = Bbend

  # - first bend starting from gun
  TEL2_bendsol3_col_zstart     = TEL2_bendsol2_col_zend + 52.9*mm
  TEL2_bendsol3_col_zend       = TEL2_bendsol3_col_zstart + 90*mm
  TEL2_bendsol3_col_length     = 90*mm
  TEL2_bendsol3_col_ri         = 193*mm
  TEL2_bendsol3_col_ro         = 265*mm
  TEL2_bendsol3_col_b          = Bbend
  
   # - Col Solenoid
  TEL2_solenoid_col_zstart     = TEL2_solenoid_main_zend + 548.26*mm
  TEL2_solenoid_col_length     = 345*mm
  TEL2_solenoid_col_zend       = TEL2_solenoid_col_zstart+TEL2_solenoid_col_length
  TEL2_solenoid_col_radi       = 120*mm
  TEL2_solenoid_col_rado       = 248*mm
  TEL2_solenoid_col_b          = Bcoll

# --- Beam size & position --- #
beama0              = 17.5e0*mm                 # Beam size in X [m]
beamb0              = 17.5e0*mm                 # Beam size in Y [m]
beamap0             = .0e0*mm                   # Beam divergence in X [m_x/m_z]
beambp0             = .0e0*mm                   # Beam divergence in Y [m_x/m_z]
beamx0              = .0e0*mm                   # Beam centroid in X [m]
beamy0              = .0e0*mm                   # Beam centroid in Y [m]
beamxp0             = .0e0*mm                   # Beam centroid velocity in X [m/s]
beamyp0             = .0e0*mm                   # Beam centroid velocity in Y [m/s]

# --- Beam inject parameters ---#
beamxinject         = .0e0*mm                   # Injected beam centroid in X [m]
beamyinject         = .0e0*mm                   # Injected beam centroid in Y [m]
beamxpinject        = .0e0*mm                   # Injected beam centroid velocity in X [m/s]
beamypinject        = .0e0*mm                   # Injected beam centroid velocity in Y [m/s]
beamainject         = 17.5*mm                   # Injected beam radius in X [m]
beambinject         = 17.5*mm                   # Injected beam radius in Y [m]
beamapinject        = .0e0*mm                   # Injected beam divergence in X [m]
beambpinject        = .0e0*mm                   # Injected beam divergence in Y [m]
beamainjmin         = 6.75*mm                   # Injected beam inner radius in X [m]
beambinjmin         = 6.75*mm                   # Injected beam inner radius in Y [m]
beamzinject         = machine_zstart

###################
# >>>  Script <<< #
###################

#------------------------------#
#     Invoke setup routine     #
#------------------------------#

setup(makepsfile=0)
winon()
palette("ImageJ_Fire.gp")

top.diposet = false


#---------------------------#
#      Particle Loading     #
#---------------------------#

# We only need to load the particles if the inject type is profle
if (machine_injtype == "profile"):
  print("Reading particle positions...")
  posi = fromfile(fns[item], sep=' ')
  npart = len(posi)/2
  posi = reshape(posi, (npart,2))
  print("Calculating charge density according to particle distribution...")

print("Number of macroparticles = %e" % npart)

#-----------------------------#
#     Particle Properties     #
#-----------------------------#

# --- Particle parameters --- #
electron_Iz             = -Current                     # current [Amps]
cyc_freq                = echarge*Bmain/emass
timestep                = compact_factor*pi/(2*cyc_freq) #50*pi/(2*cyc_freq)
electron_vz             = .0e0                      # Velocity of Electron Beam at emitting surface [m/s]
electron_ekin           = -Cathode_Potential
electron_q              = -1.e0                     # Charge state of electrons []
vthz                    = .0e0                      # Thermal Velocity of particles [m/s]
lrelativity             = true                      # Should relativisitc effects be considered [boolean]
relativity              = true                      # Level of relativistic correctness (1: scale transverse field by 1/gamma**2)

sw=int((-electron_Iz*timestep/echarge)/npart)
elec = Species(type=Electron,color=red,weight=sw)
#prot = Species(type=Proton,color=green)

elec.ibeam        = electron_Iz       # current (Amps)
print ("Beam Current: %g" % elec.ibeam)
elec.zion         = electron_q        # charge state
print ("Particle charge: %g" % elec.zion)
top.dt           = timestep          # time step size (seconds)
print ("Cyclotron Frequency: %g" % cyc_freq)
print ("Timestep: %g" % top.dt)
elec.vbeam        = electron_vz       # velocity will be calculated from ekin
print ("Particle velocity: %g" % elec.vbeam)
elec.ekin         = electron_ekin     # kinetic energy - z directed (volts)
elec.aion         = top.emass/top.amu # electrons
top.derivqty()                       # calculates vbeam from ekin
elec.lrelativ     = lrelativity
elec.relativity   = relativity

elec.vthz         = vthz
#top.vthz        = .5e0*top.vbeam*top.emit/sqrt(top.a0*top.b0) # Vthz ~ Vthperp



#ebeam=(-Cathode_voltage)*echarge+emass*clight**2
#vbeam=clight*numpy.sqrt(1-(emass*clight**2)**2/ebeam**2)
nsteps = 2.2*machine_syslen/elec.vbeam/timestep  #40    ### --- Should be the same a above. Calculate vbeam from the energy and put it in.

print("The number of time steps is: %f steps \n" % nsteps)

#---------------------#
#     Beam Design     #
#---------------------#

# - size
elec.a0 = Cathode_rado                # Beam size in X
elec.b0 = Cathode_rado               # Beam size in Y
elec.ap0 = beamap0              # Beam divergance in X
elec.bp0 = beambp0              # Beam divergance in Y
# - centroid
elec.x0 = beamx0                # Initial beam centroid in x
elec.xp0 = beamxp0              # Initial beam centroid in vx/vz
elec.y0 = beamy0                # Initial beam centroid in y
elec.yp0 = beamyp0              # Initial beam centroid in vy/vz

#-------------------#
#     Injection     #
#-------------------#

# --- Beam Injection --- #
elec.npmax = npart          # This is the maximum number of particles
top.inject = machine_emittype   # What is the limitation on the current of the beam?    
top.zinject[0] = beamzinject # location of injection source (m)

if (machine_injtype=="gun"):
  elec.npinject = int(npart**2*sw*elec.sq/(elec.ibeam*timestep*nsteps))      # number of particles to inject each time step  - NEEDS TO BE CORRECTED
  print("number of particles injected per time step: %g" % elec.npinject)
  #top.npinject = int(top.ibeam*top.dt/top.echarge)  # number of particles to inject each time step

  top.xinject[0] = beamxinject
  top.yinject[0] = beamyinject
  top.xpinject[0] = beamxpinject
  top.ypinject[0] = beamypinject
  top.ainject[0] = elec.a0
  top.binject[0] = elec.b0
  top.ainjmin[0] = Cathode_radi
  top.binjmin[0] = Cathode_radi
  top.apinject[0] = beamapinject
  top.bpinject[0] = beambpinject
  top.vzinject[0,0] = 0.0
  top.vinject[0] = -500.0       # Voltage of injector source

# --- Profile Injection --- #

if machine_injtype == 'profile':
  # >>>  add particles according to measured current density
  xinit = posi[:,0] * mm #+ Xoff*mm
  yinit = posi[:,1] * mm #+ Yoff*mm
  zinit = zeros(npart)
  vxinit = zeros(npart)
  vyinit = zeros(npart)
  vzinit = zeros(npart) + elec.vbeam
  #elec.addpart(xinit,yinit,zinit,vxinit,vyinit,vzinit)
  def hollow_cathode_source():
    if w3d.inj_js == elec.jslist[0]:
      w3d.npgrp = npart
      gchange('Setpwork3d')
      w3d.xt[:] = xinit
      w3d.yt[:] = yinit
      w3d.zt[:] = top.zinject
      w3d.uxt[:] = vthz
      w3d.uyt[:] = vthz
      w3d.uzt[:] = elec.vbeam
      # w3d.uzt[:] = top.vbeam

  installuserparticlesinjection(hollow_cathode_source)

#-----------------#
#     Lattice     #
#-----------------#

top.diposet = false

# The zero point is at the cathode 
if machine_type == "tbench":
  # - Gun Solenoid
  addnewsolenoid(zi=tbench_solenoid_gun_zstart, zf=tbench_solenoid_gun_zend, ri=tbench_solenoid_gun_radi, ro=tbench_solenoid_gun_rado, maxbz=tbench_solenoid_gun_b)
  # - Drift before main solenoid
  addnewdrft(zs=0.37, ze=0.60, ap=machine_piperad)
  # - Main Solenoid
  addnewsolenoid(zi=tbench_solenoid_main_zstart, zf=tbench_solenoid_main_zend, ri=tbench_solenoid_main_radi, ro=tbench_solenoid_main_rado, maxbz=tbench_solenoid_main_b)
  # - Drift after main solenoid
  addnewdrft(zs=2.52, ze=2.67, ap=machine_piperad)
  # - Collector Solenoid
  addnewsolenoid(zi=tbench_solenoid_col_zstart, zf=tbench_solenoid_col_zend, ri=tbench_solenoid_col_radi, ro=tbench_solenoid_col_rado, maxbz=tbench_solenoid_col_b)

elif machine_type == "TEL2s":
  # - Gun Solenoid
  addnewsolenoid(zi=TEL2_solenoid_gun_zstart, zf= TEL2_solenoid_gun_zend, ri=TEL2_solenoid_gun_radi, ro=TEL2_solenoid_gun_rado, maxbz=TEL2_solenoid_gun_b)
  # - Bend
  #addnewbend(zs=TEL2_solenoid_gun_zend,ze=TEL2_solenoid_main_zstart,rc=(TEL2_solenoid_main_zstart-TEL2_solenoid_gun_zend)/1.02) #Bend angle is 
  # - 3 Bends before main solenoid
  addnewsolenoid(zi=TEL2_bendsol1_gun_zstart, zf=TEL2_bendsol1_gun_zend, ri=TEL2_bendsol1_gun_ri, ro=TEL2_bendsol1_gun_ro, maxbz=TEL2_bendsol1_gun_b)
  addnewsolenoid(zi=TEL2_bendsol2_gun_zstart, zf=TEL2_bendsol2_gun_zend, ri=TEL2_bendsol2_gun_ri, ro=TEL2_bendsol2_gun_ro, maxbz=TEL2_bendsol2_gun_b)
  addnewsolenoid(zi=TEL2_bendsol3_gun_zstart, zf=TEL2_bendsol3_gun_zend, ri=TEL2_bendsol3_gun_ri, ro=TEL2_bendsol3_gun_ro, maxbz=TEL2_bendsol3_gun_b)
  # - Main Solenoid
  addnewsolenoid(zi=TEL2_solenoid_main_zstart, zf= TEL2_solenoid_main_zend, ri=TEL2_solenoid_main_radi, ro=TEL2_solenoid_main_rado, maxbz=TEL2_solenoid_main_b)
  # - Bend
  #addnewbend(zs=TEL2_solenoid_main_zend, ze=TEL2_solenoid_col_zstart,rc=(TEL2_solenoid_col_zstart-TEL2_solenoid_main_zend)/1.02)
  # - 3 Bends after main solenoid
  addnewsolenoid(zi=TEL2_bendsol1_col_zstart, zf=TEL2_bendsol1_col_zend, ri=TEL2_bendsol1_col_ri, ro=TEL2_bendsol1_col_ro, maxbz=TEL2_bendsol1_col_b)
  addnewsolenoid(zi=TEL2_bendsol2_col_zstart, zf=TEL2_bendsol2_col_zend, ri=TEL2_bendsol2_col_ri, ro=TEL2_bendsol2_col_ro, maxbz=TEL2_bendsol2_col_b)
  addnewsolenoid(zi=TEL2_bendsol3_col_zstart, zf=TEL2_bendsol3_col_zend, ri=TEL2_bendsol3_col_ri, ro=TEL2_bendsol3_col_ro, maxbz=TEL2_bendsol3_col_b)
  # - Collector Solenoid
  addnewsolenoid(zi=TEL2_solenoid_col_zstart, zf= TEL2_solenoid_col_zend, ri=TEL2_solenoid_col_radi, ro=TEL2_solenoid_col_rado, maxbz=TEL2_solenoid_col_b)

elif machine_type == "TEL2":
  # - Gun Solenoid
  addnewsolenoid(zi=TEL2_solenoid_gun_zstart, zf= TEL2_solenoid_gun_zend, ri=TEL2_solenoid_gun_radi, ro=TEL2_solenoid_gun_rado, maxbz=TEL2_solenoid_gun_b)
  # - Bend
  addnewbend(zs=TEL2_solenoid_gun_zend,ze=TEL2_solenoid_main_zstart,rc=(TEL2_solenoid_main_zstart-TEL2_solenoid_gun_zend)/1.02) #Bend angle is 
  # - 3 Bends before main solenoid
  addnewsolenoid(zi=TEL2_bendsol1_gun_zstart, zf=TEL2_bendsol1_gun_zend, ri=TEL2_bendsol1_gun_ri, ro=TEL2_bendsol1_gun_ro, maxbz=TEL2_bendsol1_gun_b)
  addnewsolenoid(zi=TEL2_bendsol2_gun_zstart, zf=TEL2_bendsol2_gun_zend, ri=TEL2_bendsol2_gun_ri, ro=TEL2_bendsol2_gun_ro, maxbz=TEL2_bendsol2_gun_b)
  addnewsolenoid(zi=TEL2_bendsol3_gun_zstart, zf=TEL2_bendsol3_gun_zend, ri=TEL2_bendsol3_gun_ri, ro=TEL2_bendsol3_gun_ro, maxbz=TEL2_bendsol3_gun_b)
  # - Main Solenoid
  addnewsolenoid(zi=TEL2_solenoid_main_zstart, zf= TEL2_solenoid_main_zend, ri=TEL2_solenoid_main_radi, ro=TEL2_solenoid_main_rado, maxbz=TEL2_solenoid_main_b)
  # - Bend
  addnewbend(zs=TEL2_solenoid_main_zend, ze=TEL2_solenoid_col_zstart,rc=(TEL2_solenoid_col_zstart-TEL2_solenoid_main_zend)/1.02)
  # - 3 Bends after main solenoid
  addnewsolenoid(zi=TEL2_bendsol1_col_zstart, zf=TEL2_bendsol1_col_zend, ri=TEL2_bendsol1_col_ri, ro=TEL2_bendsol1_col_ro, maxbz=TEL2_bendsol1_col_b)
  addnewsolenoid(zi=TEL2_bendsol2_col_zstart, zf=TEL2_bendsol2_col_zend, ri=TEL2_bendsol2_col_ri, ro=TEL2_bendsol2_col_ro, maxbz=TEL2_bendsol2_col_b)
  addnewsolenoid(zi=TEL2_bendsol3_col_zstart, zf=TEL2_bendsol3_col_zend, ri=TEL2_bendsol3_col_ri, ro=TEL2_bendsol3_col_ro, maxbz=TEL2_bendsol3_col_b)
  # - Collector Solenoid
  addnewsolenoid(zi=TEL2_solenoid_col_zstart, zf= TEL2_solenoid_col_zend, ri=TEL2_solenoid_col_radi, ro=TEL2_solenoid_col_rado, maxbz=TEL2_solenoid_col_b)


# >>>  Set input parameters describing the 3d simulation.
w3d.nx = 32                    # number of grid points in x
w3d.ny = 32                    # number of grid points in y
w3d.nz = 256                   # number of grid points in z  (Shouldn't we use more in z direction?)
top.prwall = machine_piperad    # radius at which particles are scraped  - 
                                  # Optional radius of cylindrical wall that absorbs particles
                                  # When zero, uses largest cylinder that fits in grid

# >>>  Set to finite beam.
w3d.xmmin = -machine_piperad    # mesh minimum in x (meters)
w3d.xmmax =  machine_piperad    # mesh maximum in x (meters)
w3d.ymmin = -machine_piperad    # mesh minimum in y (meters)
w3d.ymmax =  machine_piperad    # mesh maximum in y (meters)
w3d.zmmin =  machine_zstart      # mesh minimum in z (meters)
w3d.zmmax =  machine_syslen      # mesh maximum in z (meters)

dx = (w3d.xmmax-w3d.xmmin) / w3d.nx # Length of cell in x direction
dy = (w3d.ymmax-w3d.ymmin) / w3d.ny # Length of cell in y direction
dz = (w3d.zmmax-w3d.zmmin) / w3d.nz # Length of call in z direction

# >>>  Set up some diagnostic windows.
top.xwindows[:,1] = [-5.e-2,5.e-2]
#top.rwindows[:,1] = [0.e0,.01e0]
top.zwindows[:,1] = [machine_zstart,2*elec.vbeam*top.dt]
top.zwindows[:,2] = [machine_syslen/2-elec.vbeam*top.dt, machine_syslen/2+elec.vbeam*top.dt]
top.zwindows[:,3] = [machine_syslen-2*elec.vbeam*top.dt, machine_syslen]

# >>>  Time histories
elec.nhist = int(nsteps/10)                                             # save history data every time step
top.ifzmmnt = 2                                                 # Linear grid point weighting moments
top.itmomnts[0:3]=[0,nsteps,elec.nhist]
top.zmmntmin = machine_zstart
top.zmmntmax = machine_syslen
top.nzmmnt = w3d.nz

# # --- Setup Plots
# top.itplps[0:8]     = [0,nsteps,nsteps/10,1,5,10,20,30]
# # top.itplfreq[0:3]   = [0,nsteps,nsteps/20]
# top.itplalways[0:5] = [0,nsteps,nsteps/20,10,30]
# top.itplseldom[0:5] = [0,nsteps,nsteps/10,10,30]

#--- Setup Plots
top.itplps[0:3]     = [0,nsteps,nsteps]
# top.itplfreq[0:3]   = [0,nsteps,nsteps/20]
top.itplalways[0:3] = [0,nsteps,nsteps]
top.itplseldom[0:3] = [0,nsteps,nsteps]


#top.ipzx[0]     = seldom                                      # turn on z versus x plot
#top.ipzy[0]     = seldom
#top.iptrace[1]  = seldom                                     # turn on set of four transverse phace-space plots
#top.iptrace[2]  = seldom                                     # turn on set of four transverse phace-space plots
#top.iptrace[3]  = seldom                                     # turn on set of four transverse phace-space plots

#top.lhlinechg = false # DEFAULT:true turned of for speed until paralleldump works 
#top.lhvzofz = false # DEFAULT:true
#top.lhcurrz = false # DEFAULT:false

top.pboundnz = absorb
top.pbound0  = absorb
top.pboundxy = absorb

# >>>  set up field solver
w3d.solvergeom  = w3d.XYZgeom # use axisymmetric solver
w3d.bound0      = 1  # left boundary in z uses neumann boundary conditions
w3d.boundnz     = 1  # right boundary in z uses neumann boundary conditions
w3d.boundxy     = 0  # transverse uses dirichlet boundary conditions
if w3d.solvergeom == w3d.XYZgeom:
  # >>>  Set some flags only needed if using the 3d solver
  w3d.l4symtry    = true # turn on four-fold symmetry
  top.fstype      = 7  # use the 3d multigrid solver
  f3d.mgparam     = 1.2 # relaxation parameter for the solver
  f3d.downpasses  = 1
  f3d.uppasses    = 1
  f3d.gridmode    = 1
  f3d.mgverbose   = 1    # Default
  f3d.mgntverbose = 1  # Default
  f3d.lcndbndy    = true # Default - Turns on sub-grid boundaries
  f3d.lprecalccoeffs = true #finite difference coefficients are precalculated. Is faster but uses more memory.
  f3d.laddconductor  = false # DEFAULT Python function calladdconductor is called before fieldsolve 


top.lgridqnt = true
top.lvinject = true #false # DEFAULT, sets whether source is included in field solve. Try different setting.

# Setup Envelope Boundaries
env.zl = w3d.zmmin
env.zu = w3d.zmmax
env.dzenv = machine_syslen/1000

top.diposet = false

w3d.interpdk[1]=1
w3d.igradb=1
from loadgradb import setbsqgrad
setbsqgrad(w3d.nx,w3d.ny,w3d.nz,w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,w3d.zmmin,w3d.zmmax)





# Generate Envelope function
package("env");generate();step


# >>>  Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.).
package("w3d");generate()

# >>> Plot the PIC grid
plotgrid()

#-------------------------------#
#     Installing Conductors     #
#-------------------------------#

if (machine_injtype=="gun"):
  # --- Electron Gun --- #

  # - Gun Drift Pipe
  gun_driftpipe=ZCylinderOut(radius=Gun_pipe_ri,zlower=Gun_pipe_zstart,zupper=Gun_pipe_zend,xcent=0.0,ycent=0.0, voltage=Gun_pipe_voltage)
  # - Cathode
  gun_cathode_r = [Cathode_radi,Cathode_radi,Cathode_rado,Cathode_rado]
  gun_cathode_z = [Cathode_zstart,Cathode_zend,Cathode_zend,Cathode_zstart]
  gun_cathode_radi=[None,Cathode_radcurvb,None,None]
  gun_cathode=ZSrfrv(rsrf=gun_cathode_r,zsrf=gun_cathode_z,rad=gun_cathode_radi,voltage=Cathode_voltage,xcent=.0e0,ycent=.0e0,zcent=.0e0)
  # - Anode
  gun_anode_r           = [Anode_r1,Anode_r3,Anode_r4,Anode_r5,Anode_rendi,Anode_rendo,Anode_rendo,Anode_r2,Anode_r2,Anode_radtipo,Anode_radtipi]
  gun_anode_z           = [Anode_z1,Anode_z3,Anode_z4,Anode_z5,Anode_zend,Anode_zend,Anode_z5,Anode_z4,Anode_z2,Anode_zstart,Anode_zstart]
  gun_anode_radi        = [None,None,None,None,None,None,None,None,Anode_radcurvb,None,Anode_radcurvs]
  gun_anode             = ZSrfrv(rsrf=gun_anode_r,zsrf=gun_anode_z,rad=gun_anode_radi, voltage=Anode_voltage,xcent=.0e0,ycent=.0e0,zcent=.0e0)
  # - Electrode F
  gun_electrodef_r      = [ElectrodeF_r1,ElectrodeF_ro,ElectrodeF_ro,ElectrodeF_ri,ElectrodeF_ri]
  gun_electrodef_z      = [ElectrodeF_zend,ElectrodeF_z2,ElectrodeF_zstart,ElectrodeF_zstart,ElectrodeF_z1]
  gun_electrodef_radi   = [ElectrodeF_radcurvb,None,None,None,ElectrodeF_radcurvs]
  gun_electrodef        = ZSrfrv(rsrf=gun_electrodef_r,zsrf=gun_electrodef_z,rad=gun_electrodef_radi,voltage=ElectrodeF_voltage,xcent=0,ycent=0,zcent=.0e0)
  # - Electrode C
  gun_electrodeC_r      = [ElectrodeC_ro,ElectrodeC_ro,ElectrodeC_ri,ElectrodeC_ri]
  gun_electrodeC_z      = [ElectrodeC_z1, ElectrodeC_zstart, ElectrodeC_zstart,ElectrodeC_z1]
  gun_electrodeC_radi   = [None,None,None,ElectrodeC_radcurv]
  gun_electrodeC        = ZSrfrv(rsrf=gun_electrodeC_r,zsrf=gun_electrodeC_z,rad=gun_electrodeC_radi,voltage=ElectrodeC_voltage,xcent=0,ycent=0,zcent=.0e0)

  gun_conductors=[gun_driftpipe,gun_driftpipe,gun_cathode,gun_anode,gun_electrodef,gun_electrodeC]
  installconductor(gun_conductors)

# --- Lattice --- #

pipe = ZCylinderOut(radius=machine_piperad,zlower=Gun_pipe_zstart,zupper=machine_syslen,voltage=0.,xcent=0,ycent=0,zcent=0)
lattice_conductors=[pipe]
installconductor(lattice_conductors)

fieldsolve()  #fieldsol(-1)

#--------------------------#
#     Plotting Lattice     #
#--------------------------#

# --- Plotting Envelope Function
penv(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1,ascale=None,bscale=None,zscale=None)
fma()

# --- Plotting Potential
if (machine_injtype=="gun"):
  gun_driftpipe.draw(filled=190,color="fg")
  gun_cathode.draw(filled=160,color='fg')
  gun_anode.draw(filled=100,color='fg')
  gun_electrodef.draw(filled=150,color='fg')
  gun_electrodeC.draw(filled=150,color='fg')
pfzr(fullplane=1,plotsg=1,cond=1,fill=1,plotphi=1,plotrho=0,plotselfe=0,comp='z',titles=1)
limits(Cathode_zstart,Gun_pipe_zend,-1.2*machine_piperad,1.2*machine_piperad)
fma()

# --- Plotting Potential
if (machine_injtype=="gun"):
  gun_driftpipe.draw(filled=190,color="fg")
  gun_cathode.draw(filled=160,color='fg')
  gun_anode.draw(filled=100,color='fg')
  gun_electrodef.draw(filled=150,color='fg')
  gun_electrodeC.draw(filled=150,color='fg')
pfzr(fullplane=1,plotsg=1,cond=1,fill=1,plotphi=0,plotrho=1,plotselfe=0,comp='E',titles=1)
limits(Cathode_zstart,Gun_pipe_zend,-1.2*machine_piperad,1.2*machine_piperad)
fma()

# --- Plotting Potential
if (machine_injtype=="gun"):
  gun_driftpipe.draw(filled=190,color="fg")
  gun_cathode.draw(filled=160,color='fg')
  gun_anode.draw(filled=100,color='fg')
  gun_electrodef.draw(filled=150,color='fg')
  gun_electrodeC.draw(filled=150,color='fg')
pfzr(fullplane=1,plotsg=1,cond=1,fill=1,plotphi=0,plotrho=0,plotselfe=1,comp='E',titles=1)
limits(Cathode_zstart,Gun_pipe_zend,-1.2*machine_piperad,1.2*machine_piperad)
fma()

# --- Plotting Electric Field
if (machine_injtype=="gun"):
  gun_driftpipe.draw(filled=190,color="fg")
  gun_cathode.draw(filled=160,color='fg')
  gun_anode.draw(filled=100,color='fg')
  gun_electrodef.draw(filled=150,color='fg')
  gun_electrodeC.draw(filled=150,color='fg')
pfzr(fullplane=1,plotsg=1,cond=1,fill=1,plotphi=0,plotrho=0,plotselfe=1,comp='z',titles=1)
limits(Cathode_zstart,Gun_pipe_zend,-1.2*machine_piperad,1.2*machine_piperad)
fma()

## --- REPETITIVE PLOTS

def myplots():
  # --- Plotting Electric Field
  if (machine_injtype=="gun"):
    gun_driftpipe.draw(filled=190,color="fg")
    gun_cathode.draw(filled=160,color='fg')
    gun_anode.draw(filled=100,color='fg')
    gun_electrodef.draw(filled=150,color='fg')
    gun_electrodeC.draw(filled=150,color='fg')
    pipe.draw(filled=60,color='fg')
  limits(Cathode_zstart,machine_syslen,-1.2*machine_piperad,1.2*machine_piperad)
  # pfzr(fullplane=1,plotsg=1,cond=1,fill=1,plotphi=0,plotrho=0,plotselfe=1,comp='z',titles=1)
  ppzx(color="density",ncolor=30)
  fma()

  if (machine_injtype=="gun"):
    gun_driftpipe.draw(filled=190,color="fg")
    gun_cathode.draw(filled=160,color='fg')
    gun_anode.draw(filled=100,color='fg')
    gun_electrodef.draw(filled=150,color='fg')
    gun_electrodeC.draw(filled=150,color='fg')
    pipe.draw(filled=60,color='fg')
  limits(Cathode_zstart,machine_syslen,-1.2*machine_piperad,1.2*machine_piperad)
  # pfzr(fullplane=1,plotsg=1,cond=1,fill=1,plotphi=0,plotrho=0,plotselfe=1,comp='z',titles=1)
  ppzy(color="density",ncolor=30)
  fma()

  ppvzvperp(iw=1,color="density",ncolor=30)
  fma()

  ppvxvy(iw=1,color="density",ncolor=30)
  fma()  

  ppyvy(iw=1,color="density",ncolor=30)
  fma()  

  ppxvx(iw=1,color="density",ncolor=30)
  fma()  

  pptrace(filled=1,particles=0,contours=30)
  fma()

installplalways(myplots)

#pipe.draw(filled=150,color="fg")
#limits(Cathode_zstart,machine_syslen,-1.2*machine_piperad,1.2*machine_piperad)
#hcp()

#visualizeconductors(condid=20)
#hcp()

# >>>  run time steps and dump the results.
step(nsteps)
# step(300)

ex=getex()
ey=getey()
ez=getez()
# er=geter()
bx=getbx()
by=getby()
bz=getbz()
ke=getke()
#Find vector potential


ff = open('../Results/'+date+'/'+machine_type+"_"+date+time+"_"+machine_injtype+"_"+file_ending+'_particlepos.txt','w')
ff.write('#Warp simulation of '+machine_type+'\n# Author: Vince Moens \n# Particle Dump giving position, velocity kinetic energy and the fields at the particles. \n# Date: '+date+'\n# Time: '+time+'\n\n# Time: %10.5e s\n# Number of timesteps: %10.5e \n# Timestep size: %10.5e s\n# Solenoid Fields: %10.5e-%10.5e-%10.5e T\n# Cathode-Anode voltage: %10.5eV\n# Beam current: %10.5e A\n# Beam velocity: %10.5e m/s\n# Bem velocity over c: %10.5e \n# Kinetic Energy: %10.5e eV\n\n' %(top.dt*int(nsteps),int(nsteps),top.dt,Bgun,Bmain,Bcoll,Cathode_Potential,top.ibeam,top.vbeam,top.vbeamoc,top.ekin))
ff.write('# Number of Macroparticles: %10.5e\n# Macroparticle weight: %10.5e electrons\n# Macroparticle charge: %10.5e coulombs\n\n' %(elec.nps,elec.sw,elec.sq*elec.sw))
ff.write('X[m] Y[m] Z[m] Xv[m/s] Yv[m/s] Zv[m/s] KE[eV] E_x[V/m] E_y[V/m] E_z[V/m] B_x[T] B_y[T] B_z[T] B[T]\n')
for x,y,z,u,v,g,a,b,c,d,e,f in zip(elec.xp,elec.yp,elec.zp,elec.uxp,elec.uyp,elec.uzp,ex,ey,ez,bx,by,bz):
  ff.write('%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n' %(x,y,z,u,v,g,a,b,c,d,e,f))
ff.close()

# Obtaining electric fields on grid
allocateselfeforfieldsolve()
nx,ny,nz = array(w3d.phi.shape) #- 1
getselfe3d(w3d.phi,w3d.nxlocal,w3d.nylocal,w3d.nzlocal,
           w3d.nxguardphi,w3d.nyguardphi,w3d.nzguardphi,
           w3d.selfe,w3d.nxguarde,w3d.nyguarde,w3d.nzguarde,
           w3d.dx,w3d.dy,w3d.dz,true)
selfe = w3d.selfe[:,w3d.nxguarde:-w3d.nxguarde or None,
                    w3d.nyguarde:-w3d.nyguarde or None,
                    w3d.nzguarde:-w3d.nzguarde or None]

Ex = selfe[0,...]
Ey = selfe[1,...]
Ez = selfe[2,...]


#Provide fields for whole grid information!

ff = open('../Results/'+date+'/'+machine_type+"_"+date+time+"_"+machine_injtype+"_"+file_ending+'_Efields.txt','w')
ff.write('#Warp simulation of '+machine_type+'\n# Author: Vince Moens \n# Particle Dump \n# Date: '+date+'\n# Time: '+time+'\n\n# Time: %10.5e s\n# Number of timesteps: %10.5e \n# Timestep size: %10.5e s\n# Solenoid Fields: %10.5e-%10.5e-%10.5e T\n# Cathode-Anode voltage: %10.5eV\n# Beam current: %10.5e A\n# Beam velocity: %10.5e m/s\n# Bem velocity over c: %10.5e \n# Kinetic Energy: %10.5e eV\n\n' %(top.dt*int(nsteps),int(nsteps),top.dt,Bgun,Bmain,Bcoll,Cathode_Potential,top.ibeam,top.vbeam,top.vbeamoc,top.ekin))
ff.write('# Number of Macroparticles: %10.5e\n# Macroparticle weight: %10.5e electrons\n# Macroparticle charge: %10.5e coulombs\n\n' %(elec.nps,elec.sw,elec.sq*elec.sw))
ff.write('#Grid size in x: %10.5e\n# Grid size in y: %10.5e\n# Grid size in z: %10.5e\n# Cell size in x: %10.5e m\n# Cell size in y: %10.5e m\n# Cell size in z: %10.5e\n' %(w3d.nx,w3d.ny,w3d.nz,dx,dy,dz))
ff.write('\n\nE_x[V/m]:\n\n')
data=Ex
# Write the array to disk
# I'm writing a header here just for the sake of readability
# Any line starting with "#" will be ignored by numpy.loadtxt
ff.write('# Array shape: {0}\n'.format(data.shape))

# Iterating through a ndimensional array produces slices along
# the last axis. This is equivalent to data[i,:,:] in this case
for data_slice in data:

  # The formatting string indicates that I'm writing out
  # the values in left-justified columns 7 characters in width
  # with 2 decimal places.  
  np.savetxt(ff, data_slice, fmt='%-7.2f')

  # Writing out a break to indicate different slices...
  ff.write('# New slice\n')

ff.write('\n\nE_y[V/m]:\n\n')
data=Ey
# Write the array to disk
# I'm writing a header here just for the sake of readability
# Any line starting with "#" will be ignored by numpy.loadtxt
ff.write('# Array shape: {0}\n'.format(data.shape))

# Iterating through a ndimensional array produces slices along
# the last axis. This is equivalent to data[i,:,:] in this case
for data_slice in data:

  # The formatting string indicates that I'm writing out
  # the values in left-justified columns 7 characters in width
  # with 2 decimal places.  
  np.savetxt(ff, data_slice, fmt='%-7.2f')

  # Writing out a break to indicate different slices...
  ff.write('# New slice\n')
ff.write('\n\nE_z[V/m]:\n\n')
data=Ez
# Write the array to disk
# I'm writing a header here just for the sake of readability
# Any line starting with "#" will be ignored by numpy.loadtxt
ff.write('# Array shape: {0}\n'.format(data.shape))

# Iterating through a ndimensional array produces slices along
# the last axis. This is equivalent to data[i,:,:] in this case
for data_slice in data:

  # The formatting string indicates that I'm writing out
  # the values in left-justified columns 7 characters in width
  # with 2 decimal places.  
  np.savetxt(ff, data_slice, fmt='%-7.2f')

  # Writing out a break to indicate different slices...
  ff.write('# New slice\n')
# E_y[V/m] E_z[V/m] E[V/m] B_x[T] B_y[T] B_z[T] B[T] Ax[')
# for x,y,z,u,v,g in zip(FIND THIS):
#   ff.write('%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n' %(x,y,z,u,v,g))
ff.close()

# Obtaining magnetic fields on grid

bfield = f3d.bfield
nxguardb = bfield.nxguardb
nyguardb = bfield.nyguardb
nzguardb = bfield.nzguardb

b = bfield.b[:,nxguardb:-nxguardb or None,
                   nyguardb:-nyguardb or None,
                   nzguardb:-nzguardb or None]
Bx = b[0,...]
By = b[1,...]
Bz = b[2,...]


#Provide fields for whole grid information!

ff = open('../Results/'+date+'/'+machine_type+"_"+date+time+"_"+machine_injtype+"_"+file_ending+'_Bfields.txt','w')
ff.write('#Warp simulation of '+machine_type+'\n# Author: Vince Moens \n# Particle Dump \n# Date: '+date+'\n# Time: '+time+'\n\n# Time: %10.5e s\n# Number of timesteps: %10.5e \n# Timestep size: %10.5e s\n# Solenoid Fields: %10.5e-%10.5e-%10.5e T\n# Cathode-Anode voltage: %10.5eV\n# Beam current: %10.5e A\n# Beam velocity: %10.5e m/s\n# Bem velocity over c: %10.5e \n# Kinetic Energy: %10.5e eV\n\n' %(top.dt*int(nsteps),int(nsteps),top.dt,Bgun,Bmain,Bcoll,Cathode_Potential,top.ibeam,top.vbeam,top.vbeamoc,top.ekin))
ff.write('# Number of Macroparticles: %10.5e\n# Macroparticle weight: %10.5e electrons\n# Macroparticle charge: %10.5e coulombs\n\n' %(elec.nps,elec.sw,elec.sq*elec.sw))
ff.write('#Grid size in x: %10.5e\n# Grid size in y: %10.5e\n# Grid size in z: %10.5e\n# Cell size in x: %10.5e m\n# Cell size in y: %10.5e m\n# Cell size in z: %10.5e\n' %(w3d.nx,w3d.ny,w3d.nz,dx,dy,dz))
ff.write('\n\nB_x[T]:\n\n')
data=Bx
# Write the array to disk
# I'm writing a header here just for the sake of readability
# Any line starting with "#" will be ignored by numpy.loadtxt
ff.write('# Array shape: {0}\n'.format(data.shape))

# Iterating through a ndimensional array produces slices along
# the last axis. This is equivalent to data[i,:,:] in this case
for data_slice in data:

  # The formatting string indicates that I'm writing out
  # the values in left-justified columns 7 characters in width
  # with 2 decimal places.  
  np.savetxt(ff, data_slice, fmt='%-7.2f')

  # Writing out a break to indicate different slices...
  ff.write('# New slice\n')

ff.write('\n\nB_y[T]:\n\n')
data=By
# Write the array to disk
# I'm writing a header here just for the sake of readability
# Any line starting with "#" will be ignored by numpy.loadtxt
ff.write('# Array shape: {0}\n'.format(data.shape))

# Iterating through a ndimensional array produces slices along
# the last axis. This is equivalent to data[i,:,:] in this case
for data_slice in data:

  # The formatting string indicates that I'm writing out
  # the values in left-justified columns 7 characters in width
  # with 2 decimal places.  
  np.savetxt(ff, data_slice, fmt='%-7.2f')

  # Writing out a break to indicate different slices...
  ff.write('# New slice\n')
ff.write('\n\nB_z[T]:\n\n')
data=Bz
# Write the array to disk
# I'm writing a header here just for the sake of readability
# Any line starting with "#" will be ignored by numpy.loadtxt
ff.write('# Array shape: {0}\n'.format(data.shape))

# Iterating through a ndimensional array produces slices along
# the last axis. This is equivalent to data[i,:,:] in this case
for data_slice in data:

  # The formatting string indicates that I'm writing out
  # the values in left-justified columns 7 characters in width
  # with 2 decimal places.  
  np.savetxt(ff, data_slice, fmt='%-7.2f')

  # Writing out a break to indicate different slices...
  ff.write('# New slice\n')
# E_y[V/m] E_z[V/m] E[V/m] B_x[T] B_y[T] B_z[T] B[T] Ax[')
# for x,y,z,u,v,g in zip(FIND THIS):
#   ff.write('%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n' %(x,y,z,u,v,g))
ff.close()

dump()

os.system("mv "+machine_type+"_"+date+time+"_"+machine_injtype+"_"+file_ending+"* ../Results/"+date+"/")

