setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(jsonlite)
library(ggnewscale)
library(scales)
library(reshape2)
library(ggplot2)

# data from notebook

# ALF less
target_position_ALF_00_1 <- c(103, 256, 343, 650, 671, 747)

target_position_ALF_01_1 <- c(1, 64, 136, 137, 151, 202, 216, 217, 218, 219, 
                              220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
                              230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
                              240, 241, 242, 243, 244, 245, 246, 247, 248, 253,
                              340, 373, 396, 422, 517, 568, 845, 855, 856, 857,
                              858, 859, 892, 984, 1055, 1129)

target_position_ALF_11_1 <- c(137, 254, 271, 519, 575, 576, 577, 578, 579, 580, 
                              581, 582, 583, 584, 585, 586, 587, 588, 589, 590,
                              591, 592, 593, 594, 595, 596, 597, 598, 599, 600,
                              601, 602, 603, 604, 605, 606, 607, 608, 609, 610,
                              611, 612, 613, 614, 615, 616, 617, 618, 619, 620,
                              621, 622, 623, 624, 625, 626, 627, 628, 629, 630,
                              631, 632, 633, 634, 635, 636, 637, 638, 639, 640,
                              641, 642, 643, 644, 645, 646, 647, 648, 649, 650,
                              651, 652, 653, 654, 655, 656, 657, 658, 659, 660,
                              661, 662, 663, 664, 665, 666, 667, 668, 669, 670,
                              671, 672, 673, 674, 675, 676, 677, 678, 679, 680,
                              681, 682, 683, 684, 685, 686, 687, 688, 689, 690,
                              691, 692, 693, 694, 695, 696, 697, 698, 699, 700,
                              701, 702, 703, 704, 705, 706, 707, 708, 709, 710,
                              711, 712, 713, 714, 715, 716, 717, 718, 719, 720,
                              721, 722, 723, 724, 725, 726, 727, 728, 729, 730,
                              731, 732, 733, 734, 735, 777, 903, 937, 940, 941, 
                              1129)
# ALF more
target_position_ALF_00_0 <- c(1707, 1914, 2078, 2079, 2080, 2081, 2082, 2083, 
                              2084, 2085, 2086, 2087, 2088, 2089, 2090, 2091, 
                              2092, 2093, 2094, 2095, 2096, 2097, 2098, 2099, 
                              2100, 2101, 2102, 2103, 2104, 2105, 2106, 2107, 
                              2108, 2109, 2110, 2111, 2112, 2113, 2114, 2115, 
                              2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123,
                              2124, 2125, 2126, 2127, 2128, 2129, 2130, 2131, 
                              2132, 2133, 2134, 2135, 2136, 2137, 2138, 2139, 
                              2140, 2141, 2142, 2143, 2144, 2145, 2335)

target_position_ALF_01_0 <- c(1199, 1549, 1914, 1921, 2104)

target_position_ALF_11_0 <- c(1199, 1201, 1578)

# OA less
target_position_OA_00_1 <- c(1, 103, 136, 137, 151, 202, 216, 217, 218, 219, 
                             220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
                             230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
                             240, 241, 242, 243, 244, 245, 246, 247, 248, 254,
                             256, 271, 340, 343, 396, 422, 517, 519, 650, 671,
                             731, 747, 777, 845, 855, 856, 857, 858, 859, 892,
                             903, 937, 940, 984, 1055, 1129)

target_position_OA_01_1 <- c(845)

target_position_OA_11_1 <- c()
# OA more
target_position_OA_00_0 <- c(1199, 1201, 1549, 1707, 1914, 1921, 2335)

target_position_OA_01_0 <- c(1914)

target_position_OA_11_0 <- c(1914, 2078, 2079, 2080, 2081, 2082, 2083, 2084, 
                             2085, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 
                             2093, 2094, 2095, 2096, 2097, 2098, 2099, 2100, 
                             2101, 2102, 2103, 2104, 2105, 2106, 2107, 2108, 
                             2109, 2110, 2111, 2112, 2113, 2114, 2115, 2116, 
                             2117, 2118, 2119, 2120, 2121, 2122, 2123, 2124, 
                             2125, 2126, 2127, 2128, 2129, 2130, 2131, 2132,
                             2133, 2134, 2135, 2136, 2137, 2138, 2139, 2140,
                             2141, 2142, 2143, 2144, 2145)

# OI less
target_position_OI_00_1 <- c(1, 103, 136, 137, 202, 216, 217, 218, 219, 220, 
                             221, 222, 223, 224, 225, 226, 227, 228, 229, 230,
                             231, 232, 233, 234, 235, 236, 237, 238, 239, 240,
                             241, 242, 243, 244, 245, 246, 247, 248, 253, 254,
                             256, 271, 340, 517, 650, 671, 731, 777, 892, 903, 
                             937, 940, 984, 1055, 1129)

target_position_OI_01_1 <- c(137, 151, 396, 422, 517, 845)

target_position_OI_11_1 <- c(137, 519)
# OI more
target_position_OI_00_0 <- c(1199, 1201, 1549, 1707, 1914, 2078, 2079, 2080, 
                             2081, 2082, 2083, 2084, 2085, 2086, 2087, 2088,
                             2089, 2090, 2091, 2092, 2093, 2094, 2095, 2096, 
                             2097, 2098, 2099, 2100, 2101, 2102, 2103, 2104,
                             2105, 2106, 2107, 2108, 2109, 2110, 2111, 2112,
                             2113, 2114, 2115, 2116, 2117, 2118, 2119, 2120,
                             2121, 2122, 2123, 2124, 2125, 2126, 2127, 2128,
                             2129, 2130, 2131, 2132, 2133, 2134, 2135, 2136,
                             2137, 2138, 2139, 2140, 2141, 2142, 2143, 2144,
                             2145)

target_position_OI_01_0 <- c(1914)

target_position_OI_11_0 <- c()

# SS less
target_position_SS_00_1 <- c(1, 136, 137, 216, 217, 218, 219, 220, 221, 222, 
                             223, 224, 225, 226, 227, 228, 229, 230, 231, 232,
                             233, 234, 235, 236, 237, 238, 239, 240, 241, 242,
                             243, 244, 245, 246, 247, 248, 254, 519, 777, 855,
                             856, 857, 858, 859, 903, 937, 940, 984, 1129)

target_position_SS_01_1 <- c(103, 137, 151, 202, 253, 256, 340, 396, 422, 517, 
                             671, 747, 845, 892, 1055)

target_position_SS_11_1 <- c(271, 650, 731)
# SS more
target_position_SS_00_0 <- c(1199, 1201)

target_position_SS_01_0 <- c(1549, 1707, 1914, 1921, 2335)

target_position_SS_11_0 <- c(1914, 2078, 2079, 2080, 2081, 2082, 2083, 2084, 
                             2085, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 
                             2093, 2094, 2095, 2096, 2097, 2098, 2099, 2100,
                             2101, 2102, 2103, 2104, 2105, 2106, 2107, 2108,
                             2109, 2110, 2111, 2112, 2113, 2114, 2115, 2116,
                             2117, 2118, 2119, 2120, 2121, 2122, 2123, 2124,
                             2125, 2126, 2127, 2128, 2129, 2130, 2131, 2132,
                             2133, 2134, 2135, 2136, 2137, 2138, 2139, 2140,
                             2141, 2142, 2143, 2144, 2145)

# IG less
target_position_IG_00_1 <- c(1, 103, 136, 137, 151, 202, 253, 254, 256, 340,
                             343, 396, 422, 517, 671, 731, 777, 855, 856, 857,
                             858, 859, 892, 903, 937, 940, 984, 1055, 1129)

target_position_IG_01_1 <- c(517, 747, 845)

target_position_IG_11_1 <- c(271, 519, 650)
# IG more
target_position_IG_00_0 <- c(1199, 1201, 1707, 1914)

target_position_IG_01_0 <- c(1549, 1914, 2335)

target_position_IG_11_0 <- c(1914)

# NK less
target_position_NK_00_1 <- c(1, 137, 151, 253, 254, 256, 271, 340, 422, 517, 
                             671, 731, 777, 892, 1129)

target_position_NK_01_1 <- c(137, 202, 396, 454, 517, 1055, 1073, 1097)

target_position_NK_11_1 <- c(519, 650, 903, 940, 941)
# NK more
target_position_NK_00_0 <- c(1201)

target_position_NK_01_0 <- c(1474, 1495, 1517, 1521, 1549, 1743, 2232, 2269)

target_position_NK_11_0 <- c(1199, 1914)


# CHOOSE ONE 555
modify_matrix <- function(matrix_data, target_position_list, modified_matrix) {
  for (target_position in target_position_list) {
    closest_rows <- which(matrix_data[target_position, ] < 5)
    closest_cols <- which(matrix_data[, target_position] < 5)
    closest_both <- intersect(closest_rows, closest_cols)
    modified_matrix[closest_both, closest_both] <- matrix_data[closest_both, closest_both]
  }
  return(modified_matrix)
}
# CHOOSE ONE 333
modify_matrix <- function(matrix_data, target_position_list, modified_matrix) {
  for (target_position in target_position_list) {
    closest_rows <- which(matrix_data[target_position, ] < 3)
    closest_cols <- which(matrix_data[, target_position] < 3)
    closest_both <- intersect(closest_rows, closest_cols)
    modified_matrix[closest_both, closest_both] <- matrix_data[closest_both, closest_both]
  }
  return(modified_matrix)
}

### START FROM HERE EVERUTIME

data_0 <- fromJSON("AF3/0_NBAS/fold_nbas_full_data_0.json")
matrix_data_list_0 <- data_0["pae"]
matrix_data_0 <- matrix(unlist(matrix_data_list_0), nrow = 2371, ncol = 2371)

data_1 <- fromJSON("AF3/0_NBAS/fold_nbas_full_data_1.json")
matrix_data_list_1 <- data_1["pae"]
matrix_data_1 <- matrix(unlist(matrix_data_list_1), nrow = 2371, ncol = 2371)

data_2 <- fromJSON("AF3/0_NBAS/fold_nbas_full_data_2.json")
matrix_data_list_2 <- data_2["pae"]
matrix_data_2 <- matrix(unlist(matrix_data_list_2), nrow = 2371, ncol = 2371)

data_3 <- fromJSON("AF3/0_NBAS/fold_nbas_full_data_3.json")
matrix_data_list_3 <- data_3["pae"]
matrix_data_3 <- matrix(unlist(matrix_data_list_3), nrow = 2371, ncol = 2371)

data_4 <- fromJSON("AF3/0_NBAS/fold_nbas_full_data_4.json")
matrix_data_list_4 <- data_4["pae"]
matrix_data_4 <- matrix(unlist(matrix_data_list_4), nrow = 2371, ncol = 2371)

###
# mean
matrix_data <- (matrix_data_0 + matrix_data_1 + matrix_data_2 + 
                  matrix_data_3 + matrix_data_4) / 5

###

rownames(matrix_data) <- 1:nrow(matrix_data)
colnames(matrix_data) <- 1:ncol(matrix_data)


#### OR HERE - START FROM HERE EVERUTIME 1111111111111111111111111111111111111
modified_matrix_00_1 <- matrix(30, nrow = 2371, ncol = 2371) # ps
modified_matrix_01_1 <- matrix(30, nrow = 2371, ncol = 2371) # ps
modified_matrix_11_1 <- matrix(30, nrow = 2371, ncol = 2371) # ps
modified_matrix_00_0 <- matrix(30, nrow = 2371, ncol = 2371) # ps
modified_matrix_01_0 <- matrix(30, nrow = 2371, ncol = 2371) # ps
modified_matrix_11_0 <- matrix(30, nrow = 2371, ncol = 2371) # ps

################################################### ONLY ONE AT ONCE 222222222
################################################### ALF ALF ALF
# Modify matrices for each target 
modified_matrix_00_1 <- modify_matrix(matrix_data, target_position_ALF_00_1, modified_matrix_00_1)
modified_matrix_01_1 <- modify_matrix(matrix_data, target_position_ALF_01_1, modified_matrix_01_1)
modified_matrix_11_1 <- modify_matrix(matrix_data, target_position_ALF_11_1, modified_matrix_11_1)
# Modify matrices for each target 
modified_matrix_00_0 <- modify_matrix(matrix_data, target_position_ALF_00_0, modified_matrix_00_0)
modified_matrix_01_0 <- modify_matrix(matrix_data, target_position_ALF_01_0, modified_matrix_01_0)
modified_matrix_11_0 <- modify_matrix(matrix_data, target_position_ALF_11_0, modified_matrix_11_0)
###################################################

################################################### OPTIC ATROPHY
# Modify matrices for each target
modified_matrix_00_1 <- modify_matrix(matrix_data, target_position_OA_00_1, modified_matrix_00_1)
modified_matrix_01_1 <- modify_matrix(matrix_data, target_position_OA_01_1, modified_matrix_01_1)
modified_matrix_11_1 <- modify_matrix(matrix_data, target_position_OA_11_1, modified_matrix_11_1)
# Modify matrices for each target
modified_matrix_00_0 <- modify_matrix(matrix_data, target_position_OA_00_0, modified_matrix_00_0)
modified_matrix_01_0 <- modify_matrix(matrix_data, target_position_OA_01_0, modified_matrix_01_0)
modified_matrix_11_0 <- modify_matrix(matrix_data, target_position_OA_11_0, modified_matrix_11_0)
###################################################

################################################### OSTEOGENESIS IMPERFECTA
# Modify matrices for each target
modified_matrix_00_1 <- modify_matrix(matrix_data, target_position_OI_00_1, modified_matrix_00_1)
modified_matrix_01_1 <- modify_matrix(matrix_data, target_position_OI_01_1, modified_matrix_01_1)
modified_matrix_11_1 <- modify_matrix(matrix_data, target_position_OI_11_1, modified_matrix_11_1)
# Modify matrices for each target
modified_matrix_00_0 <- modify_matrix(matrix_data, target_position_OI_00_0, modified_matrix_00_0)
modified_matrix_01_0 <- modify_matrix(matrix_data, target_position_OI_01_0, modified_matrix_01_0)
modified_matrix_11_0 <- modify_matrix(matrix_data, target_position_OI_11_0, modified_matrix_11_0)
###################################################

################################################### SHORT STATURE
# Modify matrices for each target
modified_matrix_00_1 <- modify_matrix(matrix_data, target_position_SS_00_1, modified_matrix_00_1)
modified_matrix_01_1 <- modify_matrix(matrix_data, target_position_SS_01_1, modified_matrix_01_1)
modified_matrix_11_1 <- modify_matrix(matrix_data, target_position_SS_11_1, modified_matrix_11_1)
# Modify matrices for each target
modified_matrix_00_0 <- modify_matrix(matrix_data, target_position_SS_00_0, modified_matrix_00_0)
modified_matrix_01_0 <- modify_matrix(matrix_data, target_position_SS_01_0, modified_matrix_01_0)
modified_matrix_11_0 <- modify_matrix(matrix_data, target_position_SS_11_0, modified_matrix_11_0)
###################################################

################################################### IMMUNOGLOBULINS
# Modify matrices for each target
modified_matrix_00_1 <- modify_matrix(matrix_data, target_position_IG_00_1, modified_matrix_00_1)
modified_matrix_01_1 <- modify_matrix(matrix_data, target_position_IG_01_1, modified_matrix_01_1)
modified_matrix_11_1 <- modify_matrix(matrix_data, target_position_IG_11_1, modified_matrix_11_1)
# Modify matrices for each target
modified_matrix_00_0 <- modify_matrix(matrix_data, target_position_IG_00_0, modified_matrix_00_0)
modified_matrix_01_0 <- modify_matrix(matrix_data, target_position_IG_01_0, modified_matrix_01_0)
modified_matrix_11_0 <- modify_matrix(matrix_data, target_position_IG_11_0, modified_matrix_11_0)
###################################################

################################################### NATURAL KILLERS
# Modify matrices for each target
modified_matrix_00_1 <- modify_matrix(matrix_data, target_position_NK_00_1, modified_matrix_00_1)
modified_matrix_01_1 <- modify_matrix(matrix_data, target_position_NK_01_1, modified_matrix_01_1)
modified_matrix_11_1 <- modify_matrix(matrix_data, target_position_NK_11_1, modified_matrix_11_1)
# Modify matrices for each target
modified_matrix_00_0 <- modify_matrix(matrix_data, target_position_NK_00_0, modified_matrix_00_0)
modified_matrix_01_0 <- modify_matrix(matrix_data, target_position_NK_01_0, modified_matrix_01_0)
modified_matrix_11_0 <- modify_matrix(matrix_data, target_position_NK_11_0, modified_matrix_11_0)
###################################################


# Melt matrices for plotting 3333333333333333333333333333333333333333333333333
df_data <- melt(matrix_data)
df_modified_00_1 <- melt(modified_matrix_00_1) # bs
df_modified_01_1 <- melt(modified_matrix_01_1) # bs
df_modified_11_1 <- melt(modified_matrix_11_1) # bs
df_modified_00_0 <- melt(modified_matrix_00_0) # bs
df_modified_01_0 <- melt(modified_matrix_01_0) # bs
df_modified_11_0 <- melt(modified_matrix_11_0) # bs

# CHOOSE ONE 44444444444444444444444444444444444444444444444444444444444444444
target_positionsU <- data.frame(
  pos = c(target_position_ALF_00_1, target_position_ALF_01_1, target_position_ALF_11_1),
  color = c(rep("#0099ff", length(target_position_ALF_00_1)),
            rep("#ffd700", length(target_position_ALF_01_1)),
            rep("#E4181C", length(target_position_ALF_11_1)))
)
# CHOOSE ALF
target_positionsD <- data.frame(
  pos = c(target_position_ALF_00_0, target_position_ALF_01_0, target_position_ALF_11_0),
  color = c(rep("#0099ff", length(target_position_ALF_00_0)),
            rep("#ffd700", length(target_position_ALF_01_0)),
            rep("#E4181C", length(target_position_ALF_11_0)))
)
# CHOOSE ONE OA
target_positionsU <- data.frame(
  pos = c(target_position_OA_00_1, target_position_OA_01_1),# , target_position_OA_11_1),
  color = c(rep("#0099ff", length(target_position_OA_00_1)),
            rep("#ffd700", length(target_position_OA_01_1)))# ,
            # rep("#E4181C", length(target_position_OA_11_1)))
)
# CHOOSE OA
target_positionsD <- data.frame(
  pos = c(target_position_OA_00_0, target_position_OA_01_0, target_position_OA_11_0),
  color = c(rep("#0099ff", length(target_position_OA_00_0)),
            rep("#ffd700", length(target_position_OA_01_0)),
            rep("#E4181C", length(target_position_OA_11_0)))
)

# CHOOSE ONE OI
target_positionsU <- data.frame(
  pos = c(target_position_OI_00_1, target_position_OI_01_1, target_position_OI_11_1),
  color = c(rep("#0099ff", length(target_position_OI_00_1)),
            rep("#ffd700", length(target_position_OI_01_1)),
            rep("#E4181C", length(target_position_OI_11_1)))
)
# CHOOSE OI
target_positionsD <- data.frame(
  pos = c(target_position_OI_00_0, target_position_OI_01_0),# , target_position_OI_11_0),
  color = c(rep("#0099ff", length(target_position_OI_00_0)),
            rep("#ffd700", length(target_position_OI_01_0)))# ,
            # rep("#E4181C", length(target_position_OI_11_0)))
)

# CHOOSE ONE SS
target_positionsU <- data.frame(
  pos = c(target_position_SS_00_1, target_position_SS_01_1, target_position_SS_11_1),
  color = c(rep("#0099ff", length(target_position_SS_00_1)),
            rep("#ffd700", length(target_position_SS_01_1)),
            rep("#E4181C", length(target_position_SS_11_1)))
)
# CHOOSE SS
target_positionsD <- data.frame(
  pos = c(target_position_SS_00_0, target_position_SS_01_0, target_position_SS_11_0),
  color = c(rep("#0099ff", length(target_position_SS_00_0)),
            rep("#ffd700", length(target_position_SS_01_0)),
            rep("#E4181C", length(target_position_SS_11_0)))
)

# CHOOSE ONE IG
target_positionsU <- data.frame(
  pos = c(target_position_IG_00_1, target_position_IG_01_1, target_position_IG_11_1),
  color = c(rep("#0099ff", length(target_position_IG_00_1)),
            rep("#ffd700", length(target_position_IG_01_1)),
            rep("#E4181C", length(target_position_IG_11_1)))
)
# CHOOSE IG
target_positionsD <- data.frame(
  pos = c(target_position_IG_00_0, target_position_IG_01_0, target_position_IG_11_0),
  color = c(rep("#0099ff", length(target_position_IG_00_0)),
            rep("#ffd700", length(target_position_IG_01_0)),
            rep("#E4181C", length(target_position_IG_11_0)))
)

# CHOOSE ONE NK
target_positionsU <- data.frame(
  pos = c(target_position_NK_00_1, target_position_NK_01_1, target_position_NK_11_1),
  color = c(rep("#0099ff", length(target_position_NK_00_1)),
            rep("#ffd700", length(target_position_NK_01_1)),
            rep("#E4181C", length(target_position_NK_11_1)))
)
# CHOOSE NK
target_positionsD <- data.frame(
  pos = c(target_position_NK_00_0, target_position_NK_01_0, target_position_NK_11_0),
  color = c(rep("#0099ff", length(target_position_NK_00_0)),
            rep("#ffd700", length(target_position_NK_01_0)),
            rep("#E4181C", length(target_position_NK_11_0)))
)

# Plot
combined_plot <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.7) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_00_1, aes(fill = value)) +
  scale_fill_gradient(low = "#0099ff", high = "transparent") + # nowar
  new_scale_fill() +
  geom_tile(data = df_modified_00_0, aes(fill = value)) +
  scale_fill_gradient(low = "#0099ff", high = "transparent") + # nowar
  new_scale_fill() +
  geom_tile(data = df_modified_01_1, aes(fill = value)) +
  scale_fill_gradient(low = "#ffd700", high = "transparent") + # nowar
  new_scale_fill() +
  geom_tile(data = df_modified_01_0, aes(fill = value)) +
  scale_fill_gradient(low = "#ffd700", high = "transparent") + # nowar
  new_scale_fill() +
  geom_tile(data = df_modified_11_1, aes(fill = value)) +
  scale_fill_gradient(low = "#E4181C", high = "transparent") + # nowar, OA cancelles this
  new_scale_fill() +
  geom_tile(data = df_modified_11_0, aes(fill = value)) +
  scale_fill_gradient(low = "#E4181C", high = "transparent") + # nowar, OI cancelles this
  new_scale_fill() +
  # scale_x_continuous(name = "", breaks = seq(min(df_data$Var1), max(df_data$Var1), by = 500), labels = seq(min(df_data$Var1), max(df_data$Var1), by = 500)) +
  # scale_y_continuous(name = "", breaks = seq(min(df_data$Var2), max(df_data$Var2), by = 250), labels = seq(min(df_data$Var2), max(df_data$Var2), by = 250)) +
  theme_classic()  +
  theme(axis.text.x = element_blank(), # element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), # element_line(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

# Add vertical lines with adjustable lengths using geom_segment
combined_plot <- combined_plot +
  geom_segment(data = target_positionsD, aes(x = 0, xend = pos, y = pos, yend = pos, color = color),
               linetype = "solid", linewidth = 0.5, alpha = 0.5, show.legend = FALSE) +
  scale_color_identity()

# Add vertical lines with adjustable lengths using geom_segment
combined_plot <- combined_plot +
  geom_segment(data = target_positionsU, aes(x = 2371, xend = pos, y = pos, yend = pos, color = color),
               linetype = "solid", linewidth = 0.5, alpha = 0.5, show.legend = FALSE) +
  scale_color_identity()

# ALF
ggsave(file.path("PAE_pheno", paste("PAE_ALF.png", sep = "")), 
       combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)
# OA
ggsave(file.path("PAE_pheno", paste("PAE_OA.png", sep = "")), 
       combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)
# OI
ggsave(file.path("PAE_pheno", paste("PAE_OI.png", sep = "")), 
       combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)
# SS
ggsave(file.path("PAE_pheno", paste("PAE_SS.png", sep = "")), 
       combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)
# IG
ggsave(file.path("PAE_pheno", paste("PAE_IG.png", sep = "")), 
       combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)
# NK
ggsave(file.path("PAE_pheno", paste("PAE_NK.png", sep = "")), 
       combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)



##############################################################################
# Plot
combined_plotb <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.5) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_00, aes(fill = value)) +
  scale_fill_gradient(low = "#0099ff", high = "transparent") + # nowar
  #new_scale_fill() +
  #geom_tile(data = df_modified_01, aes(fill = value)) +
  #scale_fill_gradient(low = "#ffd700", high = "transparent") + # nowar
  #new_scale_fill() +
  #geom_tile(data = df_modified_11, aes(fill = value)) +
  #scale_fill_gradient(low = "#E4181C", high = "transparent") + # nowar
  new_scale_fill() +
  # scale_x_continuous(name = "", breaks = seq(min(df_data$Var1), max(df_data$Var1), by = 500), labels = seq(min(df_data$Var1), max(df_data$Var1), by = 500)) +
  # scale_y_continuous(name = "", breaks = seq(min(df_data$Var2), max(df_data$Var2), by = 250), labels = seq(min(df_data$Var2), max(df_data$Var2), by = 250)) +
  theme_classic()  +
  theme(axis.text.x = element_blank(), # element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), # element_line(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

# Plot
combined_plotc <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.5) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  #new_scale_fill() +
  #geom_tile(data = df_modified_00, aes(fill = value)) +
  #scale_fill_gradient(low = "#0099ff", high = "transparent") + # nowar
  new_scale_fill() +
  geom_tile(data = df_modified_01, aes(fill = value)) +
  scale_fill_gradient(low = "#ffd700", high = "transparent") + # nowar
  #new_scale_fill() +
  #geom_tile(data = df_modified_11, aes(fill = value)) +
  #scale_fill_gradient(low = "#E4181C", high = "transparent") + # nowar
  new_scale_fill() +
  # scale_x_continuous(name = "", breaks = seq(min(df_data$Var1), max(df_data$Var1), by = 500), labels = seq(min(df_data$Var1), max(df_data$Var1), by = 500)) +
  # scale_y_continuous(name = "", breaks = seq(min(df_data$Var2), max(df_data$Var2), by = 250), labels = seq(min(df_data$Var2), max(df_data$Var2), by = 250)) +
  theme_classic()  +
  theme(axis.text.x = element_blank(), # element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), # element_line(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

# Plot
combined_plotp <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.5) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  #new_scale_fill() +
  #geom_tile(data = df_modified_00, aes(fill = value)) +
  #scale_fill_gradient(low = "#0099ff", high = "transparent") + # nowar
  #new_scale_fill() +
  #geom_tile(data = df_modified_01, aes(fill = value)) +
  #scale_fill_gradient(low = "#ffd700", high = "transparent") + # nowar
  new_scale_fill() +
  geom_tile(data = df_modified_11, aes(fill = value)) +
  scale_fill_gradient(low = "#E4181C", high = "transparent") + # nowar
  new_scale_fill() +
  # scale_x_continuous(name = "", breaks = seq(min(df_data$Var1), max(df_data$Var1), by = 500), labels = seq(min(df_data$Var1), max(df_data$Var1), by = 500)) +
  # scale_y_continuous(name = "", breaks = seq(min(df_data$Var2), max(df_data$Var2), by = 250), labels = seq(min(df_data$Var2), max(df_data$Var2), by = 250)) +
  theme_classic()  +
  theme(axis.text.x = element_blank(), # element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), # element_line(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

# SS
ggsave(file.path("PAE_pheno", paste("PAE_SStestb.png", sep = "")), 
       combined_plotb,
       width = 1000, height = 1000, units = "px", dpi = 300)

# SS
ggsave(file.path("PAE_pheno", paste("PAE_SStestc.png", sep = "")), 
       combined_plotc,
       width = 1000, height = 1000, units = "px", dpi = 300)

# SS
ggsave(file.path("PAE_pheno", paste("PAE_SStestp.png", sep = "")), 
       combined_plotp,
       width = 1000, height = 1000, units = "px", dpi = 300)
##############################################################################