#' Basic models of population growth
#' 
#' Inspired by Ch. 7 of "Quantitative Ecology: A New Unified Approach" by
#'   Lehman, Loberg, and Clark.

library(shiny)
library(ggplot2)
library(patchwork)

sim_data <- function(r = 3,
                     s = -(r + 1),
                     N0 = 0.011107,
                     dt = 1,
                     t_max = 20) {
    out <- as.data.frame(matrix(NA, nrow = 0, ncol = 3))
    names(out) <- c("t", "N", "dN")
    t <- 0
    N <- N0
    while (t <= t_max) {
        dN <- N * (r + s * N) * dt
        out <- rbind(out,
                     list(t = t, N = N, dN = dN))
        # update time and population size
        t <- t + dt
        N <- N + dN
        if (N < 0) N <- 0
    }
    if (s > 0 & r > 0) {
        sing_t <- (1/r) * log(1 + r/(s * N0)) # time of singularity
    } else {
        sing_t <- Inf # no singularity
    }
    out <- within(out, {
        N <- ifelse(t <= sing_t, N, NA)
        dN_dt <- dN / dt
        growth <- dN_dt / N
    })
    
    attr(out, "r") <- r
    attr(out, "s") <- s
    attr(out, "N0") <- N0
    attr(out, "sing_t") <- sing_t
    row.names(out) <- NULL
    return(out)
}

solve_diff_eq <- function(simdat, t_max, length.out = 1000) {
    r <- attr(simdat, "r")
    s <- attr(simdat, "s")
    N0 <- attr(simdat, "N0")
    differ <- data.frame(x = seq(0, t_max, length.out = length.out))
    differ$y <- vapply(differ$x,
                       function(t) {
                           1 / (
                               ((s/r) + (1/N0)) * exp(-r * t) - (s/r)
                               )},
                          numeric(1))
    return(differ)
}

#' plot direction field
#' f: differential equation function with signature `function(t, y)`
plot_directionfield <- function(f,
                                t_range = c(-5, 5),
                                y_range = c(-5, 5),
                                radius = 0.1,
                                grid.by.t = 0.25, grid.by.y = 0.25,
                                alpha = 1){
    # credit:
    # https://stackoverflow.com/questions/47984874/how-to-create-a-slope-field-in-r
    
    # initial plot - ensure large enough
    plot(t_range, y_range,
         main = "Direction field", ylab = "y", xlab = "t",
         pch = ".")
    # plot arrows
    tlin = seq(min(t_range), max(t_range), grid.by.t)
    ylin = seq(min(y_range), max(y_range), grid.by.y)
    for(x in tlin) {
        for(y in ylin) {
            slope = f(x, y)
            if(is.na(slope)) {
                col = rgb(0, 0, 0, alpha)
            } else if(slope > 0) {
                col = rgb(0, 0, 1, alpha)
            } else if (slope < 0) {
                col = rgb(1, 0, 0, alpha)
            } else if(slope == 0) {
                col = rgb(0, 1, 0, alpha)
            }
            arrows(radius * cos(atan(slope) + pi) + x,
                   radius * sin(atan(slope) + pi) + y,
                   radius * cos(atan(slope)) + x,
                   radius * sin(atan(slope)) + y, 
                   length = 0.2 * radius, col = col)
        }
    }
}



ui <- fluidPage(
    withMathJax(),
    h1("Basic models of population growth"),
    p("Inspired by 'Quantitative Ecology: A New Unified Approach' by Lehman, Loberg, and Clark."),
    h2("Single-Species Dynamics"),
    p(paste(
        "Let N(t) describe the population size of some species evolving over time.\n",
        "The (per-capita) population growth can be described by the following differential equation,",
        "where r and s are parameters that depend on the scenario being modeled."
    )),
    uiOutput("growth_formula1"),
    tableOutput("growth_descr_table"),
    sidebarLayout(
        sidebarPanel(
            selectInput("x_axis", "X-axis",
                        choices = c("t", "N", "dN", "dN/dt", "growth per-capita"),
                        selected = "t"),
            selectInput("y_axis", "Y-axis",
                        choices = c("t", "N", "dN", "dN/dt", "growth per-capita"),
                        selected = "N"),
            sliderInput("r", "r: Intrinsic growth rate",
                        min = -5, max = 5, value = 3, step = 0.1),
            p("The intrinsic growth rate is the rate at which the population grows ",
              "in the absence of any other members of its population. It could be thought of ",
              "as the growth rate of the population after controlling for the effect of ",
              "population density, N."),
            sliderInput("s", "s: Density-dependent growth rate",
                        min = -5, max = 5, value = -4, step = 0.1),
            p("The 's' term encodes the effect of population density on the growth rate. ",
              "If s > 0, the population grows orthologistically; if s < 0, the population ",
              "grows logistically."),
            sliderInput("N0", "N0: Initial population size",
                        min = 0, max = 1, value = 0.011107, step = 0.0001),
            p("Depending on the value of other parameters, the initial population size ",
              "may have a significant effect on the population dynamics. The emergence ",
              "of very different patterns from slightly different starting points is known as ",
              "sensitive dependence on initial conditions."),
            sliderInput("dt", "dt: Time step",
                        min = 0.1, max = 1, value = 1, step = 0.1),
            p("Decreasing dt towards 0 intends to simulate infinitely small time steps, ",
              "moving from a difference equation to a differential equation."),
            sliderInput("t_max", "End period",
                        min = 10, max = 100, step = 10, value = 20)
        ),
        mainPanel(
            plotOutput("single_species_plot"),
            plotOutput("direction_field"),
            p("If s = 0, the population grows exponentially - i.e., the population is independent of its own density (N)."),
            p("If s > 0, the population grows orthologistically - i.e., the growth rate increases with population density, or the population is 'density-enhanced'. This model can only be useful at certain ranges, for it will reach a singularity point."),
            p("If s < 0, the population grows logistically - i.e., the growth rate decreases with population density, or the population is 'density-limited' and will reach a carrying capacity K = -r/s.")
        )
    )
)


server <- function(input, output, session) {
    output$growth_formula1 <- renderUI({
        withMathJax(
            "$$\\frac{1}{N} \\frac{d N}{d t} = r + s  N$$"
        )
    })
    simmed <- reactive({
        sim_data(r = as.numeric(input$r),
                 s = as.numeric(input$s),
                 N0 = as.numeric(input$N0),
                 dt = as.numeric(input$dt),
                 t_max = as.integer(input$t_max))
    })
    axes_map <- function(key) {
        c(
            "t" = "t",
            "N" = "N",
            "dN" = "dN",
            "dN/dt" = "dN_dt",
            "growth per-capita" = "growth"
        )[key]
    }
    
    subtitle_map <- function(s, r) {
        if (s > 0) {
            if (r > 0) {
                subt <- "Orthologistic growth"
            } else if (r < 0) {
                subt <- "Orthologistic decay"
            } else {
                # r = 0, then dN/dt = sN^2
                subt <- "Orthologistic growth"
            }
        } else if (s == 0) {
            if (r > 0) {
                subt <- "Exponential growth"
            } else if (r < 0) {
                subt <- "Exponential decay"
            } else { # r == 0
                subt <- "No growth"
            }
        } else { # s < 0
            if (r > 0) {
                subt <- "Logistic growth"
            } else if (r < 0) {
                subt <- "Logistic decay"
            } else {
                # r = 0, then dN/dt = -sN^2
                subt <- "Logistic decay"
            }
        }
        paste0(subt, " (r = ", r, ", s = ", s, ")")
    }
    output$single_species_plot <- renderPlot({
        x_ax <- axes_map(input$x_axis)
        y_ax <- axes_map(input$y_axis)
        dat <- simmed()
        r <- attr(dat, "r")
        s <- attr(dat, "s")
        N0 <- attr(dat, "N0")
        sing_t <- attr(dat, "sing_t")
        
        p <- ggplot(dat, aes(x = !!sym(x_ax), y = !!sym(y_ax))) +
            geom_line(alpha = 0.5) +
            geom_point(size = 1) +
            scale_color_manual(values = c("Differential Eq. Solution" = "red")) +
            labs(title = "Population Projection",
                 subtitle = subtitle_map(s, r),
                 x = input$x_axis, y = input$y_axis, color = "") +
            theme_minimal() +
            theme(legend.position = "bottom")
        
        differ <- solve_diff_eq(dat, t_max = input$t_max)
        differ[differ$x > sing_t, "y"] <- NA
        differ <- differ[!is.na(differ$y), ]
        if (x_ax == "t" & y_ax == "N") {
            p <- p +
                geom_line(data = differ,
                          aes(x = x, y = y, color = "Differential Eq. Solution"),
                          alpha = 0.5)
        }
        if (s < 0 && r > 0) {
            K <- -r/s # carrying capacity
            p <- p + labs(caption = paste("Carrying capacity K =", round(K, 3)))
        }
        if (s > 0 && r > 0) {
            p <- p +
                geom_vline(xintercept = sing_t,
                           alpha = 0.5) +
                labs(caption = paste("Singularity at t =", round(sing_t, 3)))
        }
        
        # Density Plot
        p2 <- ggplot(dat, aes(x = !!sym(y_ax))) +
            geom_histogram(bins = input$t_max/2, fill = "black") +
            labs(x = input$y_axis, y = "Count", title = "Distribution") +
            theme_minimal() +
            coord_flip(xlim = range(dat[[y_ax]]))
        
        p <- p + p2 + plot_layout(widths = c(2, 1))
        return(p)
    })
    output$growth_descr_table <- renderTable({
        tbl <- c(
            "r > 0, s > 0" = "Orthologistic growth",
            "r < 0, s > 0" = "Orthologistic growth with Allee point",
            "r > 0, s = 0" = "Exponential growth",
            "r > 0, s < 0" = "Logistic growth with a carrying capacity",
            "r < 0, s <= 0" = "Inviable population declining to extinction"
        )
        data.frame(
            "Model" = names(tbl),
            "Description" = tbl
        )
    })
    output$direction_field <- renderPlot({
        dat <- simmed()
        r <- attr(dat, "r")
        s <- attr(dat, "s")
        N0 <- attr(dat, "N0")
        y_max <- min(max(dat$N, na.rm = TRUE) + 0.1, 1)
        t_max <- input$t_max
        
        
        
        f <- \(t, N) N * (r + s * N)
        plot_directionfield(f,
                            t_range = c(0, t_max),
                            y_range = c(0, y_max),
                            grid.by.t = t_max/40,
                            grid.by.y = y_max/50,
                            radius = t_max/(200 * (t_max/20)),
                            alpha = 0.7)
    })
}


shinyApp(ui, server)

