from math import pow

print "------ Given Optimal Params Get Other Params -----"
numCellCharesPerCore = 70.0;
numComputeCharesPerCore = 14 * numCellCharesPerCore;
h = 0.04
x_chareDim = 4 * h
y_chareDim = 4 * h
z_chareDim = 4 * h
x_dim = 2.0;
y_dim = 2.0;
z_dim = 2.0;
volume = x_dim * y_dim * z_dim
chare_volume = x_chareDim * y_chareDim * z_chareDim
x_numChare = x_dim / x_chareDim
y_numChare = y_dim / y_chareDim
z_numChare = z_dim / z_chareDim

numCellChares = x_numChare * y_numChare * z_numChare
numComputeChares = 14 * numCellChares
print("Number of Cell Chares: %.10f" % numCellChares)
print("Number of Compute Chares: %.10f" % numComputeChares)

numCores = numCellChares / numCellCharesPerCore
print("Number of cores: %.10f" % numCores)


print "------ Given dimensions get optimal params -----"
# num_charesPerCore = 70.0
# num_particlesPerChare = 80.0;
# max_numCores = 120.0;
# given_r_interaction = .02;
# x_dim = 2.0;
# y_dim = 2.0;
# z_dim = 2.0;


# v_domain = x_dim * y_dim * z_dim;
# wall_padding = given_r_interaction/2;
# total_numChares = num_charesPerCore * max_numCores;
# total_numParticles = total_numChares * num_particlesPerChare;
# total_cellChares = total_numChares/3


# v_chare = v_domain / total_cellChares

# x_chareDim = pow(v_chare, 1.0/3.0);
# y_chareDim = pow(v_chare, 1.0/3.0);
# z_chareDim = pow(v_chare, 1.0/3.0);

# print("xChareDime: %.10f" % x_chareDim)
# print("yChareDime: %.10f" % y_chareDim)
# print("zChareDime: %.10f" % z_chareDim)


# calculated_num_particlesPerChareDim = (x_chareDim - (2.0 * wall_padding)) / given_r_interaction;
# given_num_particlesPerChareDim = pow(num_particlesPerChare, 1.0/3.0)

# print("Total num chares: %.10f" % total_numChares)
# print("Total num cell chares: %.10f" % total_cellChares)
# print("Total num compute chares: %.10f" % (14*total_cellChares))
# print("Volume per chare: %.10f" % v_chare)
# print("x, y, z length per chare: %.10f %.10f %.10f" % (x_chareDim, y_chareDim, z_chareDim))
# print("Calculated Num particles per chare dim [%.10f] for r = [%.10f]" % (calculated_num_particlesPerChareDim, given_r_interaction))

# calculated_r_interaction = (x_chareDim)/pow(num_particlesPerChare,1.0/3.0);

# print("Calculated interaction radius [%.10f] for [%.10f] num particles per chare" % (calculated_r_interaction, given_num_particlesPerChareDim))





