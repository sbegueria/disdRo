# Generate internal objects with disdrometer diameter and velocity bins.

# Thies _______

# particle size and velocity bin limits
dia_t <- c(0.125, 0.25, 0.375, 0.5, 0.750, 1, 1.250, 1.5, 1.75,
                      2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5)
vel_t <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.4, 1.8, 2.2, 2.6, 3, 3.4,
                      4.2, 5, 5.8, 6.6, 7.4, 8.2, 9, 10, 11)

# particle size and velocity bin widths
dia_w_t <- -(dia_t[1:22] - dia_t[2:23])
vel_w_t <- -(vel_t[1:20] - vel_t[2:21])

# particle size and velocity bins means
dia_m_t <- (dia_t[1:22] + dia_t[2:23]) / 2
vel_m_t <- (vel_t[1:20] + vel_t[2:21]) / 2



# Parsivel _______

# particle size and velocity bin limits
dia_p <- c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1,
                           1.125, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 3, 3.5, 4,
                           4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 23, 26)
vel_p <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                           1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0,
                           4.8, 5.6, 6.4, 7.2, 8.0, 9.6, 11.2, 12.8, 14.4, 16.0,
                           19.2, 22.4)

# particle size and velocity bin widths
dia_w_p <- -(dia_p[1:32] - dia_p[2:33])
vel_w_p <- -(vel_p[1:32] - vel_p[2:33])

# particle size and velocity means
dia_m_p <- (dia_p[1:32] + dia_p[2:33])/2
vel_m_p <- (vel_p[1:32] + vel_p[2:33])/2
  

# Export to R/sysdata.rda

devtools::use_data(dia_t, vel_t, dia_w_t ,vel_w_t, dia_m_t ,vel_m_t,
                   dia_p, vel_p, dia_w_p ,vel_w_p, dia_m_p ,vel_m_p,
                   internal=TRUE)
