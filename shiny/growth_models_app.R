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
        dN <- (r + s * N) * N * dt
        out <- rbind(out,
                     list(t = t, N = N, dN = dN))
        # update time and population size
        t <- t + dt
        N <- N + dN
    }
    out <- within(out, {
        dN_dt <- dN / dt
        growth <- dN_dt / N
    })
    attr(out, "r") <- r
    attr(out, "s") <- s
    attr(out, "N0") <- N0
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


ui <- fluidPage(
    withMathJax(),
    h1("Basic models of population growth"),
    p("Inspired by 'Quantitative Ecology: A New Unified Approach' by Lehman, Loberg, and Clark."),
    p("H. Elliott"),
    h2("Single-Species Dynamics"),
    uiOutput("growth_formula1"),
    sidebarLayout(
        sidebarPanel(
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
                        min = 20, max = 100, step = 20, value = 20),
            selectInput("x_axis", "X-axis",
                        choices = c("t", "N", "dN", "dN/dt", "growth"),
                        selected = "t"),
            selectInput("y_axis", "Y-axis",
                        choices = c("t", "N", "dN", "dN/dt", "growth"),
                        selected = "N")
        ),
        mainPanel(
            plotOutput("single_species_plot"),
            tableOutput("growth_descr_table"),
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
    output$single_species_plot <- renderPlot({
        x_ax <- input$x_axis
        y_ax <- input$y_axis
        if (x_ax == "dN/dt") x_ax <- "dN_dt"
        if (y_ax == "dN/dt") y_ax <- "dN_dt"
        
        p <- ggplot(simmed(), aes_string(x = x_ax, y = y_ax)) +
            geom_line(alpha = 0.5) +
            geom_point(size = 1) +
            scale_color_manual(values = c("Differential Eq. Solution" = "red")) +
            labs(title = "Population Projection",
                 x = x_ax, y = y_ax, color = "") +
            theme_minimal() +
            theme(legend.position = "bottom")
        
        differ <- solve_diff_eq(simmed(), t_max = input$t_max)
        if (x_ax == "t" & y_ax == "N") {
            p <- p +
                geom_line(data = differ,
                          aes(x = x, y = y, color = "Differential Eq. Solution"),
                          alpha = 0.5)
        }
        
        r <- attr(simmed(), "r")
        s <- attr(simmed(), "s")
        N0 <- attr(simmed(), "N0")
        # if (s < 0 && r > 0) {
        #     K <- -r/s # carrying capacity
        #     p <- p + geom_hline(yintercept = K, linetype = "dashed")
        # }
        # if (s > 0) {
        #     t <- (1/r) * log(1 + r/(s * N0)) # time of singularity
        #     p <- p + geom_vline(xintercept = t, linetype = "dashed")
        # }
        
        # Density Plot
        p2 <- ggplot(subset(simmed(), t > 5),
                     aes_string(x = "N")) +
            geom_histogram(bins = 10) +
            labs(x = "N", y = "Count", title = "Distribution after t = 5") +
            theme_minimal() +
            coord_flip()
        
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
}


shinyApp(ui, server)

