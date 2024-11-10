
#' Quantile function, it has as purpose to estimate the 
#' quantile given an symbolic_histogram 
#' @param x it is a symbolic_histogram object
#' @param prob it is the probability of the quantile to estimate
#' @return return a numeric value
quantile_symbolic_histogram <- function(x, prob = 0.1){
#    print(x)
    breaks <- x$breaks
    props <- x$props 

    if(prob != 1){
        
        w_i <- c(0, cumsum(props))
        idx <- which(prob < w_i)[1]
        idx <- idx - 1
        w_l <- w_i[idx]
        w_r <- w_i[idx + 1]

        z_l <- breaks[idx]
        z_r <- breaks[idx + 1]

        z_l + (prob - w_l)/(w_r - w_l)*(z_r - z_l)
    }else{
        max(breaks)
    }
    

}
quantile_symbolic_histogram <- Vectorize(quantile_symbolic_histogram)








#' Function to generate uniform histograms
#' @param x symbolic_histogram object that is going to be uniformed 
#' @param cum_w cummulative weight that is considered to uniform histogram. It can be NULL, but is used length_quantile in this case
#' @param length_quantiles distance beteween the weights to uniform histogram. It can be NULL, but is used cum_w in this case
#' @return symbolic histogram uniformed
#' @keywords internal
#' @export
uniform_histogram <- function(x, cum_w = NULL,  length_quantiles = 0.1){

    if(is.null(cum_w) & is.null(length_quantiles)) stop("You must introduce cum_w, or length_quantiles")

    if(is.null(cum_w)){
        w <- seq(0, 1, length_quantiles)
        props <- rep(length_quantiles, length(w)-1)
    }else{
        w <- cum_w 
        props <- w[-1] - w[-length(w)]
        #props <- props[!is.na(props)]
    }
    

    breaks <- quantile_symbolic_histogram(x,  w)

    out <- list(
        breaks = breaks, 
        props = props
    )

    vctrs::new_vctr(list(out), class = "symbolic_histogram")
    
}





#' Function to estimate the mean of symbolic_histograms
#' @param x vector of symbolic_histogram objects
#' @return a symbolic_histogram objects that is the Frechet mean by Mallows
#' @export 
mean.symbolic_histogram <- function(x){



    if(!all(sapply(x$props, function(y) identical(y, x$props)[1]))) stop("Proportions are different, please uniform the histograms")
    if(!all(sapply(presion_diastolica$props, length) == lengths(presion_diastolica$props)[1])) stop("Different number of bins in the histograms, please uniform the histograms")


    breaks <- x$breaks 

    breaks_mean <- Reduce(`+`, breaks)/length(x)
    props <- x$props[1]
    out <- list(
        breaks = breaks_mean, 
        props = props
    )
    vctrs::new_vctr(list(out), class = "symbolic_histogram")   
}


#' Function to estimate the median of symbolic_histograms
#' @param x symbolic_histogram objects
#' @return the meadian of the symbolic_histogram
#' @export 
median.symbolic_histogram <- function(x){
    
    breaks <- x$breaks 

    props <- x$props


    w_i <- c(0, cumsum(x$props))
    idx <- which(w_i > 0.5)[1]
    idx <- idx-1

    
    q <- breaks[idx] + 0.5/(w_i[idx+1] - w_i[idx]) * (breaks[idx+1] - breaks[idx])

    q

}
median.symbolic_histogram <- Vectorize(median.symbolic_histogram)



#' This function estimate the distance between 
#' two object of kind symbolic_histogram, of package 
#' RSDA
#' @param x a symbolic_histogram object
#' @param y a symbolic_histogram object
#' @return return the Mallows Distance between the objects
#' @export
d_mallows <- function(x, y){
    breaks_x <- x$breaks
    breaks_y <- y$breaks
    length_x <- length(breaks_x)-1
    length_y <- length(breaks_y)-1

    if(length_x != length_y) stop(glue::glue("Different number of bins, x has {length_x} while y has {length_y}"))

    props_x <- x$props
    props_y <- y$props

    if(!identical(props_x, props_y)) stop(glue::glue("Different proportions, x has {props_x} while y has {props_y}"))

    centers_x <- (breaks_x[-length(breaks_x)] + breaks_x[-1])/2

    

    centers_y <- (breaks_y[-length(breaks_y)] + breaks_y[-1])/2

    ranges_x <- (-breaks_x[-length(breaks_x)] + breaks_x[-1])


    ranges_y <- (-breaks_y[-length(breaks_y)] + breaks_y[-1])

    
    # print(length(props_x))
    # print(length(centers_x-centers_y))
    sum_centers <- props_x*(centers_x-centers_y)^2

    sum_ranges <- props_x*1/3*(ranges_x-ranges_y)^2

    return(sqrt(sum(sum_centers) + sum(sum_ranges)))

}

d_mallows <- Vectorize(d_mallows)





# presion_diastolica <- c(
#     vctrs::new_vctr(list(list(breaks = c(50,60,70), props = c(0.4, 0.6))), class = "symbolic_histogram"), 
#     vctrs::new_vctr(list(list(breaks = c(70,80,90), props = c(0.2, 0.8))), class = "symbolic_histogram"),
#     vctrs::new_vctr(list(list(breaks = c(90,92,100), props = c(0.5, 0.5))), class = "symbolic_histogram"),
#     vctrs::new_vctr(list(list(breaks = c(80,85,108), props = c(0.6, 0.4))), class = "symbolic_histogram"),
#     vctrs::new_vctr(list(list(breaks = c(50,63,70), props = c(0.4, 0.6))), class = "symbolic_histogram"),
#     vctrs::new_vctr(list(list(breaks = c(80,90,100), props = c(0.5, 0.5))), class = "symbolic_histogram"),
#     vctrs::new_vctr(list(list(breaks = c(60,80,100), props = c(0.2, 0.8))), class = "symbolic_histogram"),
#     vctrs::new_vctr(list(list(breaks = c(76,85,90), props = c(0.5, 0.5))), class = "symbolic_histogram"),
#     vctrs::new_vctr(list(list(breaks = c(70,100,110), props = c(0.4, 0.6))), class = "symbolic_histogram"),
#     vctrs::new_vctr(list(list(breaks = c(90,100,110), props = c(0.4, 0.6))), class = "symbolic_histogram"),
#     vctrs::new_vctr(list(list(breaks = c(78,88,100), props = c(0.2, 0.8))), class = "symbolic_histogram")
#     )


# df_cardiological_hist <- data.frame(
#     individuo = c(1:11)
# )


# df_cardiological_hist$presion_diastolica <- presion_diastolica


