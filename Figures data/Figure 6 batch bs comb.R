setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(jsonlite)
library(ggnewscale)
library(scales)
library(reshape2)
library(ggplot2)

data <- fromJSON("AF-A2RRP1-F1-predicted_aligned_error_v4.json")

matrix_data_list <- data[["predicted_aligned_error"]]
matrix_data <- matrix(unlist(matrix_data_list), nrow = 2371, ncol = 2371)

rownames(matrix_data) <- 1:nrow(matrix_data)
colnames(matrix_data) <- 1:ncol(matrix_data)

target_position_list_RAB18 <- c(1077, 1078, 1081, 1082, 1085, 1128, 1130, 1176, 1179,
                                1180, 1181, 1182, 1183, 1184, 1185, 1186, 1190, 1192, 
                                1231, 1232, 1233, 1234, 1236, 1237)

target_position_list_USE1 <- c(727, 730, 731, 732, 733, 734, 735, 736, 737, 761, 
                               762, 763, 764, 765, 766, 767, 768, 771, 821, 822, 
                               824, 825, 826, 852, 853, 854, 855, 856, 857, 858, 
                               859, 860, 861, 862, 863, 864, 865, 888, 889, 890, 
                               891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 
                               901, 902, 903, 921, 926)

target_position_list_ZW10 <- c(1305, 1307, 1308, 1309, 1310, 1311, 1313, 1315, 
                               1316, 1319, 1320, 1323, 1334, 1335, 1336, 1337, 
                               1338, 1339, 1340, 1341, 1342, 1343, 1344, 1345, 
                               1346, 1347, 1348, 1349, 1350, 1351, 1352, 1353, 
                               1354, 1355, 1356, 1357, 1358, 1359, 1360, 1361, 
                               1362, 1364, 1365, 1464, 1466, 1467, 1468, 1469, 
                               1470, 1471, 1472, 1473, 1474, 1475, 1476, 1477, 
                               1478, 1479, 1480, 1481, 1482, 1483, 1484, 1485, 
                               1487, 1488, 1489, 1490, 1491, 1492, 1493, 1494, 
                               1495, 1496, 1497, 1498, 1499, 1500, 1501, 1502, 
                               1503, 1504, 1505, 1506, 1507, 1508, 1509, 1510, 
                               1511, 1512, 1513, 1514, 1515, 1516, 1517, 1518, 
                               1519, 1520, 1521, 1522, 1523, 1524, 1525, 1526, 
                               1527, 1528, 1529, 1530, 1532, 1533, 1535, 1536, 
                               1537, 1539, 1540, 1541, 1542, 1543, 1544, 1545, 
                               1546, 1547, 1548, 1549, 1550, 1551, 1552, 1553, 
                               1554, 1556, 1557, 1565, 1566, 1567, 1568, 1569, 
                               1570, 1571, 1572, 1573, 1574, 1575, 1576, 1577, 
                               1579, 1580, 1581, 1582, 1583, 1584, 1585, 1586, 
                               1587, 1588, 1589, 1590, 1591, 1592, 1593, 1594, 
                               1595, 1596, 1597, 1598, 1599, 1600, 1601, 1602, 
                               1603, 1604, 1605, 1606, 1607, 1608, 1609, 1610, 
                               1611, 1612, 1613, 1614, 1615, 1616, 1617, 1618, 
                               1619, 1620, 1621, 1622, 1624, 1625, 1629, 1635, 
                               1636, 1638, 1639, 1640, 1641, 1642, 1643, 1644, 
                               1645, 1646, 1647, 1648, 1649, 1650, 1651, 1657, 
                               1661, 1681, 1684, 1685, 1686, 1687, 1688, 1689, 
                               1691, 1692)

target_position_list_UPF3B <- c(1652, 1653, 1654, 1655, 1656, 1666, 1669, 1670, 
                                1672, 1673, 1674, 1676, 1703, 1706, 1707, 1710, 
                                1711, 1746)

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
modified_matrix_3 <- matrix(NA, nrow = 2371, ncol = 2371) # change to NA
modified_matrix_4 <- matrix(NA, nrow = 2371, ncol = 2371) # change to NA
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
modified_matrix_1 <- modify_matrix(matrix_data, target_position_list_RAB18, modified_matrix_1)
modified_matrix_2 <- modify_matrix(matrix_data, target_position_list_USE1, modified_matrix_2)
modified_matrix_3 <- modify_matrix(matrix_data, target_position_list_ZW10, modified_matrix_3)
modified_matrix_4 <- modify_matrix(matrix_data, target_position_list_UPF3B, modified_matrix_4)
modified_matrix_5 <- modify_matrix_pv(matrix_data, all_affected_aa, modified_matrix_5)

####################################################
# Find row and column numbers with values less than 30
rows_less_than_30 <- which(rowSums(modified_matrix_5 < 30) > 0)
cols_less_than_30 <- which(colSums(modified_matrix_5 < 30) > 0)

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
df_modified_4 <- melt(modified_matrix_4)
df_modified_5 <- melt(modified_matrix_5)

# Rescale values for each modified matrix
rescale_values <- function(df_modified) {
  df_modified$rescaled_value <- rescale(df_modified$value)
  return(df_modified)
}

# Merge the two data frames
df_modified_combined <- merge(df_modified_3, df_modified_4, by = c("Var1", "Var2"))

# Find paired NA values and fill them with 30 in both value.x and value.y columns
df_modified_combined$value.x[is.na(df_modified_combined$value.x) & !is.na(df_modified_combined$value.y)] <- 30
df_modified_combined$value.y[is.na(df_modified_combined$value.y) & !is.na(df_modified_combined$value.x)] <- 30

df_modified_combined$rescaled_value.x <- rescale(df_modified_combined$value.x)
df_modified_combined$rescaled_value.y <- rescale(df_modified_combined$value.y)

###
combined_colors <- vector("character", length = nrow(df_modified_combined))

combined_colors[!is.na(df_modified_combined$rescaled_value.x) & 
                  !is.na(df_modified_combined$rescaled_value.y)] <- rgb(
                    0, 
                    (1 - df_modified_combined$rescaled_value.x[!is.na(df_modified_combined$rescaled_value.x) & 
                                                                       !is.na(df_modified_combined$rescaled_value.y)]), 
                    (1 - df_modified_combined$rescaled_value.y[!is.na(df_modified_combined$rescaled_value.x) & 
                                                                  !is.na(df_modified_combined$rescaled_value.y)]))

# Assign NA to the rows where NA values were present
combined_colors[is.na(combined_colors)] <- NA

# Assign the combined colors to the combined value column
df_modified_combined$combined_value <- combined_colors

# Replace empty (missing) values with NA
df_modified_combined$combined_value[df_modified_combined$combined_value == ""] <- NA


# Plot
combined_plot <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.7) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_1, aes(fill = value)) +
  scale_fill_gradient(low = "purple", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_2, aes(fill = value)) +
  scale_fill_gradient(low = "yellow", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_combined, aes(fill = combined_value)) +
  scale_fill_identity(na.value = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_5, aes(fill = value)) +
  scale_fill_gradient(low = "red", high = "transparent") +
  new_scale_fill() +
  theme_classic() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")

#
# Plot PV
combined_plot <- ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(data = df_data, aes(fill = value), alpha = 0.7) +
  scale_fill_gradient(low = "darkgreen", high = "transparent") +
  new_scale_fill() +
  geom_tile(data = df_modified_5, aes(fill = value)) +
  scale_fill_gradient(low = "red", high = "transparent") +
  new_scale_fill() +
  theme_classic() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")

#
#
ggsave(file.path("PAE/test_bs", paste("BSSS.png", sep = "")), 
       combined_plot,
       width = 1000, height = 1000, units = "px", dpi = 300)

