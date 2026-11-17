#installing packages
if (!require(remotes)) {
  install.packages("remotes")
}
remotes::install_github("evo-palaeo/treesurgeon")
library(treesurgeon)

#ape - analyses of phylogenetics and evolution
install.packages("ape")
library(ape)

# visualisation tool
if (!require(ggplot2)) {
  install.packages("ggplot2")
}

library(ggplot2)

if (!require(cluster)) {
  install.packages("cluster")
}
library(cluster)

library(dplyr)

# Import aligned FASTA file
sequences <- read.dna("~/Documents/Bioinformatics MSc/Bioinformatics/Mystery Meat/sequences_by_outlet/working_directory/aligned_sequences.fasta", format = "fasta")

#Examine the data
sequences
image(
  sequences,
  show.labels = TRUE,
  cex = 0.4,
  col = c(hcl.colors(n = 4, palette = "zissou"), "grey", "white")
)

# Calculatye Pairwise Sequence Distances
raw_distances <- dist.dna(
  sequences,
  model = "raw"
  )
as.matrix(raw_distances)

# Performing a PCoA
jc_dist <- ape::dist.dna(
  sequences,
  model = "JC69"
)

pcoa_result <- cmdscale(
  jc_dist,
  k = 10,
  eig = TRUE
)

#extract PCoA coordinates
pcoa_df <- as.data.frame(pcoa_result$points)

#name the coordinate axes
colnames(pcoa_df) <- paste0("PCoA", seq_len(ncol(pcoa_df)))

#add sequence IDs
pcoa_df$Sequence_ID <- rownames(pcoa_df)

#extract food outlet names
pcoa_df$Outlet <- sub(
  "_seq.*$",
  "",
  pcoa_df$Sequence_ID)

table(pcoa_df$Outlet)

# Calculate percentages from positive eigenvalues
eig <- pcoa_result$eig
variance <- 100 * eig / sum(eig[eig > 0])

# Plot PCoA Plots

# Plot PCoA1 against PCoA2
p1 <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = Outlet)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PCoA1 (", round(variance[1], 1), "%)"),
    y = paste0("PCoA2 (", round(variance[2], 1), "%)"),
    title = "PCoA: Axes 1 and 2",
    colour = "Outlet"
  ) +
  theme_classic()

p1

p2 <- ggplot(pcoa_df, aes(x = PCoA3, y = PCoA4, colour = Outlet)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PCoA3 (", round(variance[3], 1), "%)"),
    y = paste0("PCoA4 (", round(variance[4], 1), "%)"),
    title = "PCoA: Axes 3 and 4",
    colour = "Outlet"
  ) +
  theme_classic()

p2

# Extract all 10 PCoA axes for clustering
pcoa_coords <- pcoa_df[, paste0("PCoA", 1:10)]

# Perform k-means clustering

km <- kmeans(
  pcoa_coords,
  centers = 3,
  nstart = 25
)

# Store the cluster assignments
pcoa_df$Cluster <- factor(km$cluster)

# Calculate Euclidean distances in the PCoA space
pcoa_dist <- dist(pcoa_coords)

# Test different numbers of clusters
k_values <- 2:10
mean_silhouette <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  
  km_test <- kmeans(
    pcoa_coords,
    centers = k_values[i],
    nstart = 25
  )
  
  sil <- silhouette(km_test$cluster, pcoa_dist)
  
  mean_silhouette[i] <- mean(sil[, "sil_width"])
}

# Plot the results
plot(
  k_values,
  mean_silhouette,
  type = "b",
  pch = 19,
  xlab = "Number of clusters (k)",
  ylab = "Mean silhouette width",
  main = "Choosing the Number of Clusters"
)

# Select the number of clusters with the highest mean silhouette width
optimal_k <- k_values[which.max(mean_silhouette)]

# Rerun k-means using the optimal number of clusters
set.seed(123)

km_optimal <- kmeans(
  pcoa_coords,
  centers = optimal_k,
  nstart = 25
)

# Store cluster assignments
pcoa_df$Clusters <- factor(km_optimal$cluster)

#Points colured by kmeans cluster
ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = Clusters)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    x = paste0("PCoA1 (", round(variance[1], 1), "%)"),
    y = paste0("PCoA2 (", round(variance[2], 1), "%)"),
    title = "K-means Clustering of Isolated Bacterial Samples",
    colour = "Cluster"
  ) +
  theme_classic()

# Second PCoA Plot
pcoa_df$Sequence_Group <- factor(km_optimal$cluster)

p3 <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = Sequence_Group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PCoA1 (", round(variance[1], 1), "%)"),
    y = paste0("PCoA2 (", round(variance[2], 1), "%)"),
    title = "PCoA: Inferred bacterial sequence groups",
    colour = "Sequence Group"
  ) +
  theme_classic()
p3

p4 <- ggplot(pcoa_df, aes(x = PCoA3, y = PCoA4, colour = Sequence_Group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PCoA3 (", round(variance[3], 1), "%)"),
    y = paste0("PCoA4 (", round(variance[4], 1), "%)"),
    title = "PCoA: Inferred bacterial sequence groups",
    colour = "Sequence Group"
  ) +
  theme_classic()
p4

# Compare the number of sequences per group at each outlet
table(pcoa_df$Outlet, pcoa_df$Sequence_Group)

#Visualising distribution of sequence groups across food outlets
ggplot(pcoa_df, aes(
  x = Outlet,
  fill = factor(Sequence_Group)
)) +
  geom_bar(position = "fill") +
  labs(
    title = "Distribution of Sequence Groups Across Food Outlets",
    x = "Food Outlet",
    y = "Proportion of Sequences",
    fill = "Sequence Group"
  ) +
  theme_minimal()


#Finding the middle point of each cluster from inferred sequence groups
centroids <- pcoa_df |>
  group_by(Sequence_Group) |>
  summarise(
    centroid_x = mean(PCoA1),
    centroid_y = mean(PCoA2)
  )

centroids

#Calculate the euclidean distance and then use slice_min to identify the smallest distance from the centroid
#Calculate Euclidean Distance
  #d = sqrt((x2 - x1)^2 + (y2 - y1)^2)
reps <- pcoa_df |>
  left_join(centroids, by = "Sequence_Group") |>
  mutate(
    distance = sqrt(
      (PCoA1 - centroid_x)^2 +
        (PCoA2 - centroid_y)^2
    )
  )

#Find Sequences' Coordinates with Smallest Distance
reps <- reps |>
  group_by(Sequence_Group) |>
  slice_min(distance, n = 1, with_ties = FALSE)

#What are these sequences?
reps |>
  select(Sequence_Group, Sequence_ID)
