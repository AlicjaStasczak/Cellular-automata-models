# Cellular automaton model
# ---------------------------------------------------
# This script runs a stochastic cellular automaton model using experimental
# cell-count data. The model logic follows the original
# working version: cells move, die, divide, and the simulated final population
# is accepted when it falls within the predefined deviation range from the
# experimental final cell count.
#
# Required input structure:
#   project_dir/
#   ├── exp_data/            
#   ├── reference.rda           # Fitted logspline distributions:
#   │                        # fit.angles, fit.distances, fit.cycle.length
#   └── this_script.R
#
# Main output:
#   simulation_result.csv

options(repos = c(CRAN = "https://cran.r-project.org"))

# -----------------------------------------------------------------------------
# 1. Project path
# -----------------------------------------------------------------------------
experimental_data_dir <- "exp_data/"
initial_conditions_file <- "reference.rda"
output_file_csv <- "simulations_results.csv"

# -----------------------------------------------------------------------------
# 2. Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "logspline", "openxlsx", "TeachingDemos", "ggplot2",
  "ggpubr", "DT", "waiter", "foreach", "doParallel", "bslib",
  "stringr", "readxl", "glue", "phylolm", "RColorBrewer",
  "globals", "ape"
)

installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

invisible(lapply(required_packages, library, character.only = TRUE, quietly = TRUE))

# -----------------------------------------------------------------------------
# 3. Model parameters
# -----------------------------------------------------------------------------

model <- "generations & neighbourhood"
cell_size <- 200
cell_matrix_size <- 70
simulation_time <- 164
deviation_percent <- 10
max_iterations <- 1000

# -----------------------------------------------------------------------------
# 4. Helper functions
# -----------------------------------------------------------------------------

# Load all files from a selected folder.
load_multiple_excel_files <- function(folder_path, sheet_name = "Overall") {
  files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)

  data_list <- lapply(files, function(file) {
    openxlsx::read.xlsx(file, sheet = sheet_name, skipEmptyRows = TRUE)
  })

  return(list(data_list = data_list, files = files))
}

# Euclidean distance between two vectors of the same length.
euclid_dist <- function(a, b) {
  sqrt(sum((a - b)^2))
}

# Angle between two 2D vectors.
angle_val <- function(M, N) {
  atan2(N[2], N[1]) - atan2(M[2], M[1])
}

# Calculate new cell coordinates based on distance and angle.
coord_val <- function(x, y, distance, angle) {
  data.frame(
    x = x + distance * cos(angle),
    y = y + distance * sin(angle)
  )
}

# Generate an empty cellular automaton grid.
generate_mesh <- function(cell_matrix_size) {
  matrix("0", nrow = cell_matrix_size, ncol = cell_matrix_size)
}

# Calculate a sequence of movement angles from x/y coordinates.
angles <- function(x, y) {
  angle <- c()

  for (i in 1:(length(x) - 1)) {
    dx <- x[i + 1] - x[i]
    dy <- y[i + 1] - y[i]
    angle[i] <- atan2(dy, dx) - atan2(0, 1)
  }

  return(angle)
}

# Calculate coordinates for multiple cells.
coord <- function(x, y, distance, angle) {
  coord_matrix <- matrix(nrow = length(x), ncol = 2)

  for (i in 1:length(x)) {
    coord_matrix[i, 1] <- x[i] + distance[i] * cos(angle[i])
    coord_matrix[i, 2] <- y[i] + distance[i] * sin(angle[i])
  }

  return(coord_matrix)
}

# Generate a random sample from a fitted density object.
rden <- function(n, den) {
  diffs <- diff(den$x)
  stopifnot(all(abs(diff(den$x) - mean(diff(den$x))) < 1e-9))

  total <- sum(den$y)
  den$y <- den$y / total
  ydistr <- cumsum(den$y)
  yunif <- runif(n)
  indices <- sapply(yunif, function(y) min(which(ydistr > y)))
  x <- den$x[indices]

  return(x)
}

# -----------------------------------------------------------------------------
# 5. Optional reference-model preparation
# -----------------------------------------------------------------------------
# The original workflow used pre-fitted logspline distributions stored in
# reference.rda. The commented function below documents how such a file can be
# generated from reference data, but it is not executed by default.

# data_parameters <- function(file_path, data_ref) {
#   data_raw <- data_ref
#   data_raw <- data_raw[!is.na(as.numeric(data_raw$Default_Labels)), ]
#
#   data <- openxlsx::read.xlsx(file_path, skipEmptyRows = TRUE)
#   data <- data[!is.na(as.numeric(data$Default_Labels)), ]
#
#   cycle_length <- data$Cycle_length
#   fit.cycle.length <- logspline::logspline(cycle_length)
#
#   dist <- as.numeric(data_raw$Dist_data)
#   fit.distances <- logspline::logspline(dist)
#
#   angle <- angles(x = data_raw$Cell.Position.X, y = data_raw$Cell.Position.Y)
#   fit.angles <- logspline::logspline(angle)
#
#   save(fit.angles, fit.distances, fit.cycle.length, file = "reference.rda")
# }

# -----------------------------------------------------------------------------
# 6. Cellular automaton initialization
# -----------------------------------------------------------------------------

automaton_init <- function(
    cell_size,
    cell_matrix_size,
    sim_time,
    cell_number_init,
    deviation_percent,
    data_exp,
    initial_conditions_path = initial_conditions_file
) {
  real_data <- data_exp
  real_data <- real_data[-(1:2), ]
  real_data <- head(real_data, -2)

  total_cells <- as.numeric(real_data[nrow(real_data), 2])
  initial_cells_experiment <- cell_number_init

  deviation_fraction <- deviation_percent / 100
  population_dev_down <- round(total_cells - (total_cells * deviation_fraction))
  population_dev_up <- round(total_cells + (total_cells * deviation_fraction))

  population_info <- data.frame(
    total.cells = total_cells,
    init.cells = initial_cells_experiment,
    lowest.population = population_dev_down,
    highest.population = population_dev_up
  )

  cell_matrix <- generate_mesh(cell_matrix_size)

  load(file = initial_conditions_path)
  unit <- cell_size / cell_matrix_size

  cell_data <- data.frame(
    x = NA,
    y = NA,
    cell.cycle = NA,
    live.time = NA,
    status = factor(NA, levels = c("Alive", "Dead", "Divided"))
  )

  for (cell_num in 1:cell_number_init) {
    s_x <- sample(1:cell_matrix_size, 1)
    s_y <- sample(1:cell_matrix_size, 1)
    x1 <- s_x * unit
    y1 <- s_y * unit
    cell_matrix[s_x, s_y] <- cell_num

    cell_data[as.character(cell_num), ] <- data.frame(
      x = x1,
      y = y1,
      cell.cycle = logspline::rlogspline(1, fit.cycle.length),
      live.time = 0,
      status = "Alive"
    )
  }

  return(list(cell_data = cell_data, cell_matrix = cell_matrix, population_info = population_info))
}

# -----------------------------------------------------------------------------
# 7. Cellular automaton simulation
# -----------------------------------------------------------------------------

run_cellular_automaton <- function(
    data_exp,
    input_file_name,
    s,
    cell_size,
    cell_matrix_size,
    sim_time,
    cell_number_init,
    cell_cycle_shortening,
    deviation_percent,
    max_iter = NULL,
    model,
    initial_conditions_path = initial_conditions_file,
    output_file = output_file_csv
) {
  death_probability <- 0.01

  load(file = initial_conditions_path)

  real_data <- data_exp
  real_data <- real_data[-(1:2), ]
  real_data <- head(real_data, -2)

  total_cells <- as.numeric(real_data[nrow(real_data), 2])
  initial_cells_experiment <- cell_number_init

  deviation_fraction <- deviation_percent / 100
  min0 <- round(total_cells - (total_cells * deviation_fraction))
  max0 <- round(total_cells + (total_cells * deviation_fraction))

  if (is.null(cell_number_init)) {
    cell_number_init <- initial_cells_experiment
  }

  sim_iteration <- 1
  success1 <- FALSE
  out_of_range <- FALSE
  last_iteration_results <- NULL

  while (success1 == FALSE && out_of_range == FALSE) {
    unit <- cell_size / cell_matrix_size
    cell_matrix <- generate_mesh(cell_matrix_size)

    distance_values <- logspline::rlogspline(5000, fit.distances)
    angle_values <- logspline::rlogspline(5000, fit.angles)
    cycle_time_values <- logspline::rlogspline(5000, fit.cycle.length)

    time <- 0:sim_time

    cell_data <- data.frame(
      x = NA,
      y = NA,
      cellcycle = NA,
      live.time = NA,
      status = factor(NA, levels = c("Alive", "Dead", "Divided"))
    )

    for (cell_num in 1:cell_number_init) {
      s_x <- sample(1:cell_matrix_size, 1)
      s_y <- sample(1:cell_matrix_size, 1)
      x1 <- s_x * unit
      y1 <- s_y * unit
      cell_matrix[s_x, s_y] <- cell_num

      cell_data[as.character(cell_num), ] <- data.frame(
        x = x1,
        y = y1,
        cellcycle = logspline::rlogspline(1, fit.cycle.length),
        live.time = 0,
        status = "Alive"
      )
    }

    cell_number <- numeric(sim_time + 1)
    cell_stack <- list()

    for (t in time) {
      cell_matrix_new <- matrix("0", nrow = cell_matrix_size, ncol = cell_matrix_size)
      isamp <- sample(1:cell_matrix_size)
      jsamp <- sample(1:cell_matrix_size)

      for (i in isamp) {
        for (j in jsamp) {
          if (cell_matrix[i, j] != "0") {
            dead_cell <- sample(1:(total_cells - t * death_probability), 1)
            id <- cell_matrix[i, j]

            if (dead_cell == 1) {
              cell_data[id, "status"] <- "Dead"
            }

            if (!(id %in% rownames(cell_data))) {
              cat("\n=====================================\n")
              cat("ERROR: cell is missing from cell_data\n")
              cat("ID =", id, "\n")
              cat("time =", t, " i =", i, " j =", j, "\n")
              stop("Cell does not exist in cell_data")
            }

            if (!is.na(cell_data[id, "status"]) && cell_data[id, "status"] != "Dead") {
              if (cell_data[id, 3] > cell_data[id, 4]) {
                # Cell movement and ageing.
                cell_data[id, 4] <- cell_data[id, 4] + 1

                success <- FALSE
                dens <- 0

                while (success == FALSE) {
                  dens <- dens + 1
                  x <- cell_data[id, 1]
                  y <- cell_data[id, 2]
                  distance <- logspline::rlogspline(1, fit.distances) + dens * 0.05
                  angle <- logspline::rlogspline(1, fit.angles)
                  coord <- coord_val(x, y, distance, angle)
                  coord[coord < 0] <- 0.01
                  coord[coord > cell_size] <- cell_size - 0.01
                  xnew <- ceiling(coord$x / unit)
                  ynew <- ceiling(coord$y / unit)

                  if (cell_matrix_new[xnew, ynew] == "0") {
                    cell_data[id, 1] <- coord$x
                    cell_data[id, 2] <- coord$y
                    cell_matrix_new[xnew, ynew] <- id
                    success <- TRUE
                  }
                }
              } else {
                # Cell division.
                for (cell_divided in 1:2) {
                  success <- FALSE
                  dens <- 0

                  while (success == FALSE) {
                    dens <- dens + 1
                    x <- cell_data[id, 1]
                    y <- cell_data[id, 2]
                    distance <- logspline::rlogspline(1, fit.distances) + dens * 0.05
                    angle <- logspline::rlogspline(1, fit.angles)
                    coord <- coord_val(x, y, distance, angle)
                    coord[coord < 0] <- 0.01
                    coord[coord > cell_size] <- cell_size - 0.01
                    xnew <- ceiling(coord$x / unit)
                    ynew <- ceiling(coord$y / unit)

                    if (cell_matrix_new[xnew, ynew] == "0") {
                      generation <- nchar(paste(id, as.character(cell_divided), sep = ""))

                      if (model == "generations") {
                        cellcycle <- logspline::rlogspline(1, fit.cycle.length) -
                          (generation - 1) * cell_cycle_shortening
                      }

                      if (model == "neighbourhood") {
                        cellcycle <- logspline::rlogspline(1, fit.cycle.length)
                        nh_matrix <- cell_matrix[
                          max((xnew - 2), 0):min((xnew + 2), cell_matrix_size),
                          max((ynew - 2), 0):min((ynew + 2), cell_matrix_size)
                        ]
                        r_cycle <- sum(nh_matrix != "0") / prod(dim(nh_matrix))
                        cellcycle <- cellcycle - r_cycle * cell_cycle_shortening
                      }

                      if (model == "generations & neighbourhood") {
                        g_cycle <- logspline::rlogspline(1, fit.cycle.length) -
                          (generation - 1) * cell_cycle_shortening
                        nh_matrix <- cell_matrix[
                          max((xnew - 2), 0):min((xnew + 2), cell_matrix_size),
                          max((ynew - 2), 0):min((ynew + 2), cell_matrix_size)
                        ]
                        r_cycle <- sum(nh_matrix != "0") / prod(dim(nh_matrix))
                        cellcycle <- g_cycle - r_cycle * 25 * cell_cycle_shortening
                      }

                      new_cell <- data.frame(
                        x = coord$x,
                        y = coord$y,
                        cellcycle = cellcycle,
                        live.time = 0,
                        status = "Alive"
                      )

                      new_id <- paste0(id, cell_divided)
                      rownames(new_cell) <- new_id
                      cell_data <- rbind(cell_data, new_cell)
                      cell_matrix_new[xnew, ynew] <- new_id

                      success <- TRUE
                    }
                  }

                  cell_data[id, 5] <- "Divided"
                }
              }
            }
          }
        }
      }

      cell_matrix <- cell_matrix_new
      cell_stack[[t + 1]] <- cell_matrix

      cell_number[t + 1] <- sum(cell_matrix != "0")
      message(glue::glue(
        "Cell number: {cell_number[t + 1]}, time: {t}, iteration: {sim_iteration}, ",
        "Total: {total_cells}, Minimal: {min0}, Maximal: {max0}, ",
        "file: {input_file_name}, SI: {s}"
      ))

      if (cell_number[t + 1] > 0) {
        current_result <- data.frame(
          cell_number = cell_number[t + 1],
          time = t,
          iteration = sim_iteration,
          s = s,
          output_file = input_file_name
        )
      }

      if (length(cell_number[t + 1]) > 0) {
        if (cell_number[t + 1] == 0 || cell_number[t + 1] > max0) {
          last_iteration_results <- data.frame(
            cell_number = cell_number[t + 1],
            time = t,
            iteration = sim_iteration,
            s = s,
            output_file = input_file_name
          )
          break
        }
      }
    }

    if (!is.na(cell_number[sim_time + 1])) {
      if ((cell_number[sim_time + 1] >= min0) && (cell_number[sim_time + 1] <= max0)) {
        last_iteration_results <- data.frame(
          cell_number = cell_number[t + 1],
          time = t,
          iteration = sim_iteration,
          s = s,
          output_file = input_file_name
        )
        success1 <- TRUE
      }
    }

    if (!is.null(max_iter) && sim_iteration == max_iter) {
      last_iteration_results <- data.frame(
        cell_number = cell_number[t + 1],
        time = t,
        iteration = sim_iteration,
        s = s,
        output_file = input_file_name
      )
      out_of_range <- TRUE
    }

    if (is.null(max_iter) && sim_iteration == 1000) {
      last_iteration_results <- data.frame(
        cell_number = cell_number[t + 1],
        time = t,
        iteration = sim_iteration,
        s = s,
        output_file = input_file_name
      )
      out_of_range <- TRUE
    } else {
      sim_iteration <- sim_iteration + 1
    }
  }

  if (!is.null(last_iteration_results)) {
    write.table(
      last_iteration_results,
      file = output_file,
      sep = ";",
      row.names = FALSE,
      col.names = !file.exists(output_file),
      append = file.exists(output_file)
    )
  }

  return(list(sim_iteration = sim_iteration, cell_matrix = cell_matrix))
}

# -----------------------------------------------------------------------------
# 8. Batch execution
# -----------------------------------------------------------------------------

run_algorithm <- function(s, repeats = 2) {
  loaded_data <- load_multiple_excel_files(experimental_data_dir)
  experimental_data_list <- loaded_data$data_list
  files <- loaded_data$files

  for (repeat_id in 1:repeats) {
    for (file_id in seq_along(experimental_data_list)) {
      data_exp <- experimental_data_list[[file_id]]
      input_file_name <- tools::file_path_sans_ext(basename(files[file_id]))

      real_data <- data_exp
      real_data <- real_data[-(1:2), ]
      real_data <- head(real_data, -2)
      initial_cell_number <- as.numeric(real_data[1, 2])

      run_cellular_automaton(
        data_exp = data_exp,
        input_file_name = input_file_name,
        s = s,
        cell_size = cell_size,
        cell_matrix_size = cell_matrix_size,
        sim_time = simulation_time,
        cell_number_init = initial_cell_number,
        cell_cycle_shortening = s,
        deviation_percent = deviation_percent,
        max_iter = max_iterations,
        model = model,
        initial_conditions_path = initial_conditions_file,
        output_file = output_file_csv
      )
    }
  }
}

# Run the algorithm for SI values with a step of 0.1.
for (s in seq(0, 1, by = 0.1)) {
  run_algorithm(s)
}
