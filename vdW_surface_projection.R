# INFO ----------------------------
# XYZ BOLEAN UNIT 
# Before using the code optimise structures (do not use structures directly from PubChem, ect)
# Get xyz coordinate file, the file should be in .txt format for import
# Record any ring Atoms (or atoms where the plane should be placed for Bolean unit calculation)


# Required packages ----------------------------

install.packages("sf")
install.packages("ggplot2")
library('sf')
library('ggplot2')


# IMPORT DATA ----------------------------

file_name <- file.choose() #file name stored as character too for saving file 

# SANITY CHECK (molecule structure) ----------------------------
lines <- readLines(file_name, warn = FALSE, encoding = "UTF-8")

parts <- strsplit(trimws(lines), "\\s+")

is_coordinate_line <- vapply(
  parts,
  function (x) {
    length(x) >= 4 &&
      x[1] %in% c("H","C","N") &&
      !is.na(suppressWarnings(as.numeric(x[2]))) &&
      !is.na(suppressWarnings(as.numeric(x[3]))) &&
      !is.na(suppressWarnings(as.numeric(x[4])))
  },
  logical(1)
)

coordinate_parts <- parts[is_coordinate_line]

atoms <- data.frame(
  
  element = vapply(
    coordinate_parts,
    '[',
    character(1), 1),
  
  x = as.numeric(vapply(
    coordinate_parts,
    '[',
    character(1),
    2)),
  
  y = as.numeric(vapply(
    coordinate_parts, 
    '[',
    character(1),
    3)),
  
  z = as.numeric(vapply(
    coordinate_parts,
    '[',
    character(1),
    4))
)

print(atoms)
cat("\nNumber of atoms:", nrow(atoms), "\n")
cat("Elment counts:\n")
print(table(atoms$element))

# GENERAL SETTINGS ----------------------------

#define plane_atoms

coordinates <-as.matrix(atoms[, c("x", "y", "z")])

plane_coordinates  <- coordinates[plane_atoms, ]

vdW_radii <- c(H = 1.2, C=1.7, N=1.55)

#TODO add more atoms 

#define specimen_name 



# PLANE FITTING ----------------------------
plane_centre <- colMeans(plane_coordinates)
centered_plane <- sweep(plane_coordinates, 2, plane_centre, "-")

plane_svd <- svd(centered_plane)

x_axis <- plane_svd$v[,1]
normal_vector <- plane_svd$v[,3]
y_axis <- c(
  normal_vector[2] * x_axis[3] - normal_vector[3] * x_axis[2],
  
  normal_vector[3] * x_axis[1] - normal_vector[1] * x_axis[3],
  
  normal_vector[1] * x_axis[2] - normal_vector[2] * x_axis[1]
)

y_axis <- y_axis / sqrt(sum(y_axis^2))

cat("\nPlane normal vector:\n")
print(normal_vector)

centred_coordinates <- sweep(coordinates, 2, plane_centre, "-")

projected_x <- as.vector(centred_coordinates %*% x_axis)
projected_y <- as.vector(centred_coordinates %*% y_axis)

atoms$projected_x <- projected_x
atoms$projected_y <- projected_y

print(
  atoms[, c("element", "projected_x", "projected_y")]
)


# BONDI vdW RADII ----------------------------

atoms$radius <- unname(vdW_radii[atoms$element])

print(atoms[, c("element", "projected_x", "projected_y", "radius")])

atomic_centres <- st_as_sf(
  atoms,
  coords = c("projected_x", "projected_y"),
  remove = FALSE
)

atomic_circles <- st_buffer(
  atomic_centres,
  dist = atoms$radius,
  nQuadSegs = 1000
)


# BOOLEAN Union ----------------------------

footprint_union <- st_union(atomic_circles)
union_area_A2 <- as.numeric(st_area(footprint_union))

cat("\nOverlap-corrected Boolean union area:",
    union_area_A2,
    "ring(A)^2\n")
    # per cation\n")

summed_area_A2 <- sum(pi*atoms$radius^2)
overlap_area_A2 <- summed_area_A2 - union_area_A2
overlap_percent <- (overlap_area_A2 / summed_area_A2) * 100
union_area_m2 <- union_area_A2 * 1e-20

cat("\nCalculation summary\n")
cat("---------------------\n")

cat("Uncorrected sum of atomic circles:",
    sprintf("%.6f", summed_area_A2),
    "ring(A)^2\n")

cat("Boolean union footprint:",
    sprintf("%.6f", union_area_A2),
    "ring(A)^2\n")
    #/cation\n")

cat("Overlapping area removed:",
    sprintf("%.6f", overlap_area_A2),
    "ring(A)^2\n")

cat("Percentage removed:",
    sprintf("%.2f", overlap_percent),
    "%\n")

cat("Converted footprint:",
    format(union_area_m2, scientific = TRUE, digits = 7),
    "m^2\n")
    #/cation\n")



# PROJECTIONS ----------------------------
ggplot() +
  geom_sf(
    data = atomic_circles,
    aes(fill = element),
    colour = "grey30",
    alpha = 0.2,
    linewidth = 0.2
  ) +
  geom_sf(
    data = footprint_union,
    fill = NA,
    colour = "red",
    linewidth = 1.1
  ) + 
  geom_point(
    data = atoms,
    aes(x = projected_x,
        y = projected_y,
        colour = element),
    size = 3
  ) +
  coord_sf(datum = NA) +
  scale_fill_manual(
    values = c(
      C= "grey30",
      N = "royalblue",
      H = "grey85"
    )
  ) +
  scale_colour_manual(
    values = c( C="black", N="blue", H="grey50")
  ) +
  labs(
    title = paste("Projected van der Waals footprint", specimen_name),    #*" [" * C[2] * "mim]" ^"+"),
    subtitle = paste0("Boolean union area = ", 
                      round(union_area_A2,2), 
                      expression("Å²")),
                                 #/cation")),
    x = "Projected coordinate (Å)",
    y = "Projected coordinate (Å)",
    fill = "Element",
    colour = "Element"
  ) +
  theme_classic(base_size = 20)

# SAVE (image with vdW projection) ----------------------------
#The file will be saved at the working directory
#To see working directory run:
#gtwd()

ggsave(filename =  paste(specimen_name, "_vdW_projection.png"),
       width = 7,
       height = 8,
       dpi = 600)
