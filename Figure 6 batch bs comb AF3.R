setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(jsonlite)
library(ggnewscale)
library(scales)
library(reshape2)
library(ggplot2)

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

target_position_list_RAB18 <- c(1077, 1078, 1081, 1085, 1128, 1179, 1180, 1183,
                                1184, 1186, 1231, 1233, 1234, 1237)

target_position_list_USE1 <- c(734, 764, 765, 766, 856, 857, 858, 859, 861, 892,
                               893, 899, 900)

target_position_list_ZW10 <- c(1336, 1338, 1339, 1340, 1342, 1343, 1347, 1348, 
                               1350, 1351, 1354, 1357, 1358, 1470, 1473, 1474, 
                               1475, 1476, 1477, 1478, 1479, 1480, 1481, 1490,
                               1491, 1492, 1493, 1495, 1496, 1497, 1501, 1505,
                               1508, 1509, 1512, 1515, 1516, 1542, 1543, 1546,
                               1547, 1549, 1550, 1552, 1568, 1569, 1572, 1573, 
                               1576, 1584, 1585, 1586, 1587, 1588, 1590, 1591,
                               1592, 1593, 1594, 1595, 1596, 1597, 1598, 1599,
                               1600, 1601, 1602, 1603, 1605, 1608, 1612, 1640,
                               1643)

all_affected_aa <- c(256, 1, 517, 519, 777, 1549, 270, 271, 1053, 1055, 2335,
                     803, 804, 1578, 2346, 1073, 2104, 568, 64, 1097, 842, 
                     845, 340, 343, 857, 348, 95, 1121, 103, 1129, 873, 877,
                     373, 1914, 892, 1918, 1921, 134, 903, 136, 137, 650, 
                     396, 151, 409, 1178, 414, 671, 422, 937, 426, 1707, 
                     940, 941, 1708, 1199, 1201, 2232, 187, 447, 448, 1474,
                     454, 202, 1743, 1495, 984, 731, 2269, 227, 747, 1517,
                     2287, 1521, 253, 254)

modified_matrix_1 <- matrix(30, nrow = 2371, ncol = 2371)
modified_matrix_2 <- matrix(30, nrow = 2371, ncol = 2371)
modified_matrix_3 <- matrix(30, nrow = 2371, ncol = 2371)
modified_matrix_5 <- matrix(30, nrow = 2371, ncol = 2371)

# Modify matrices for each target position list
modify_matrix <- function(matrix_data, target_position_list, modified_matrix) {
  for (target_position in target_position_list) {
    closest_rows <- which(matrix_data[target_position, ] < 5)
    closest_cols <- which(matrix_data[, target_position] < 5)
    closest_both <- intersect(closest_rows, closest_cols)
    modified_matrix[closest_both, closest_both] <- matrix_data[closest_both, closest_both]
  }
  return(modified_matrix)
}

# FILTER WUTH 3 FOR PATHOGEN VAR ----- CHOOSE
modify_matrix_pv <- function(matrix_data, target_position_list, modified_matrix) {
  for (target_position in target_position_list) {
    closest_rows <- which(matrix_data[target_position, ] < 3)
    closest_cols <- which(matrix_data[, target_position] < 3)
    closest_both <- intersect(closest_rows, closest_cols)
    modified_matrix[closest_both, closest_both] <- matrix_data[closest_both, closest_both]
  }
  return(modified_matrix)
}
# FILTER WUTH 5 FOR PATHOGEN VAR ----- CHOOSE
modify_matrix_pv <- function(matrix_data, target_position_list, modified_matrix) {
  for (target_position in target_position_list) {
    closest_rows <- which(matrix_data[target_position, ] < 5)
    closest_cols <- which(matrix_data[, target_position] < 5)
    closest_both <- intersect(closest_rows, closest_cols)
    modified_matrix[closest_both, closest_both] <- matrix_data[closest_both, closest_both]
  }
  return(modified_matrix)
}


# Modify matrices for each target position list
modified_matrix_1 <- modify_matrix(matrix_data, target_position_list_USE1, modified_matrix_1)
modified_matrix_2 <- modify_matrix(matrix_data, target_position_list_RAB18, modified_matrix_2)
modified_matrix_3 <- modify_matrix(matrix_data, target_position_list_ZW10, modified_matrix_3)
modified_matrix_5 <- modify_matrix_pv(matrix_data, all_affected_aa, modified_matrix_5)

####################################################
# Find row and column numbers with values less than 30
rows_less_than_30 <- which(rowSums(modified_matrix_3 < 30) > 0)
cols_less_than_30 <- which(colSums(modified_matrix_3 < 30) > 0)

# Combine row and column numbers into a single vector
all_numbers <- c(rows_less_than_30, cols_less_than_30)

# Get unique numbers
unique_numbers <- unique(all_numbers)

output <- paste(unique_numbers, collapse = ",")
cat(output)
# Print the unique numbers
# print(unique_numbers)
######################################################

# Melt matrices for plotting
df_data <- melt(matrix_data)
df_modified_1 <- melt(modified_matrix_1)
df_modified_2 <- melt(modified_matrix_2)
df_modified_3 <- melt(modified_matrix_3)
df_modified_5 <- melt(modified_matrix_5)

# Plot
combined_plot <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.7) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_1, aes(fill = value)) +
  scale_fill_gradient(low = "#0099ff", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_2, aes(fill = value)) +
  scale_fill_gradient(low = "purple", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_3, aes(fill = value)) +
  scale_fill_gradient(low = "#ffd700", high = "transparent") +
  # new_scale_fill() +
  # geom_tile(data = df_modified_5, aes(fill = value)) +
  # scale_fill_gradient(low = "#E4181C", high = "transparent") +
  new_scale_fill() +
  theme_classic() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")

####################################################################
target_position_list_RAB18V <- c(1186, 1231, 1233, 1234, 1237)
target_position_list_RAB18H <- c(1077, 1078, 1081, 1085, 1128, 1179, 1180, 1183,
                                1184)

# Combine all target positions into one data frame
target_positionsV <- data.frame(
  pos = c(target_position_list_RAB18V, target_position_list_ZW10),
  color = c(rep("purple", length(target_position_list_RAB18V)),
            rep("#ffd700", length(target_position_list_ZW10)))
)
# Combine all target positions into one data frame
target_positionsH <- data.frame(
  pos = c(target_position_list_USE1, target_position_list_RAB18H),
  color = c(rep("#0099ff", length(target_position_list_USE1)),
            rep("purple", length(target_position_list_RAB18H)))
)

# Create the combined plot
combined_plot <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.7) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_1, aes(fill = value)) +
  scale_fill_gradient(low = "#0099ff", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_2, aes(fill = value)) +
  scale_fill_gradient(low = "purple", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_3, aes(fill = value)) +
  scale_fill_gradient(low = "#ffd700", high = "transparent") +
  theme_classic() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")

# Add vertical and horizontal lines with adjustable lengths using geom_segment
combined_plot <- combined_plot +
  geom_segment(data = target_positionsV, aes(x = pos, xend = pos, y = 0, yend = pos, color = color),
               linetype = "solid", linewidth = 0.1, alpha = 0.5, show.legend = FALSE) +
  scale_color_identity()

# Add vertical and horizontal lines with adjustable lengths using geom_segment
combined_plot <- combined_plot +
  geom_segment(data = target_positionsH, aes(x = pos, xend = pos, y = pos, yend = 2371, color = color),
               linetype = "solid", linewidth = 0.1, alpha = 0.5, show.legend = FALSE) +
  scale_color_identity()

# Save the plot
ggsave(file.path("PAE/test_bs", "AF3-lines_BS.png"), 
       plot = combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)

#######################################################################
####################################################################

sorted_positions <- sort(all_affected_aa)
cat(paste(sorted_positions, collapse = ", "))

all_affected_aaV <- c(1199, 1201, 1474, 1495, 1517, 1521, 1549, 1578, 
                      1707, 1708, 1743, 1914, 1918, 1921, 2104, 2232, 2269, 
                      2287, 2335, 2346)

all_affected_aaH <- c(1, 64, 95, 103, 134, 136, 137, 151, 187, 202, 227, 253,
                      254, 256, 270, 271, 340, 343, 348, 373, 396, 409, 414, 
                      422, 426, 447, 448, 454, 517, 519, 568, 650, 671, 731,
                      747, 777, 803, 804, 842, 845, 857, 873, 877, 892, 903,
                      937, 940, 941, 984, 1053, 1055, 1073, 1097, 1121, 1129,
                      1178)
# Combine all target positions into one data frame
target_positionsV <- data.frame(
  pos = c(all_affected_aaV),
  color = c(rep("#E4181C", length(all_affected_aaV)))
)
# Combine all target positions into one data frame
target_positionsH <- data.frame(
  pos = c(all_affected_aaH),
  color = c(rep("#E4181C", length(all_affected_aaH)))
)
# Plot PV
combined_plot <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.7) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_5, aes(fill = value)) +
  scale_fill_gradient(low = "#E4181C", high = "transparent") +
  new_scale_fill() +
  theme_classic() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")

# Add vertical and horizontal lines with adjustable lengths using geom_segment
combined_plot <- combined_plot +
  geom_segment(data = target_positionsV, aes(x = pos, xend = pos, y = 0, yend = pos, color = color),
               linetype = "solid", linewidth = 0.1, alpha = 0.5, show.legend = FALSE) +
  scale_color_identity()

# Add vertical and horizontal lines with adjustable lengths using geom_segment
combined_plot <- combined_plot +
  geom_segment(data = target_positionsH, aes(x = pos, xend = pos, y = pos, yend = 2371, color = color),
               linetype = "solid", linewidth = 0.1, alpha = 0.5, show.legend = FALSE) +
  scale_color_identity()

# Save the plot
ggsave(file.path("PAE/test_bs", "AF3-lines_PV.png"), 
       plot = combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)

#
# Plot PV
combined_plot <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.7) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_5, aes(fill = value)) +
  scale_fill_gradient(low = "#E4181C", high = "transparent") +
  new_scale_fill() +
  theme_classic() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")

############################################################
# THIS WAS USED FOR SUBDOMAIN BORDERS CLARIFICATION 
# 475
# 640
# 1380
# 1660
# 2120
line_coord = 2120
ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.7) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  new_scale_fill() +
  geom_vline(xintercept = line_coord, color = "black", linetype = "dashed", linewidth = 0.5) +  # Use linewidth instead of size
  # geom_hline(yintercept = line_coord, color = "black", linetype = "dashed", linewidth = 0.5) +  # Use linewidth instead of size
  theme_classic() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")
############################################################
ggsave(file.path("PAE/test_bs", paste("AF3-PV5.png", sep = "")), 
       combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)

