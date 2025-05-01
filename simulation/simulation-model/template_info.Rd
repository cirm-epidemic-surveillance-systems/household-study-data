\docType{class}
\title{Augmented Onsets with Single Introduction}
\name{housepi template : Augmented Onsets with Single Introduction}
\description{
MISSING DESCRIPTION.
Missing data and asymptomatic status are not modeled.
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
\itemize{
\item Force of infection: \deqn{
{\lambda}_{i}(t) = {\mu}^{(\rm sus)}_{i} {\beta} n^{-{\gamma}} \sum_{k\,:\,t_{k}<t_{i}} {\mu}^{(\rm inf)}_{k} f_{k}(t) \, {\mu}^{(\rm con)}_{k{\to}i}
}{
lambda_i(t) = mu_sus_i * beta * n^(-gamma) sum_{k: t_k < t_i} mu_inf_k f_k(t) mu_con_ki
}
\item Time of infection: \eqn{t_{i}}{t_i} = \code{onset_date} - \code{incub_period}
\item Infectivity profile: \eqn{f_{k}}{f_k}
\item `Susceptibility` \eqn{{\mu}^{(\rm sus)}}{mu_sus}: 
\tabular{rrr}{
                            |\tab \strong{FALSE}   |\tab \strong{TRUE}      \cr
        \strong{vaccinated} |\tab \emph{reference} |\tab \code{sus_vacc}
}
\item `Infectivity` \eqn{{\mu}^{(\rm inf)}}{mu_inf}: 
\tabular{rrr}{
                            |\tab \strong{FALSE}   |\tab \strong{TRUE}       \cr
        \strong{vaccinated} |\tab \emph{reference} |\tab \code{inf_vacc} 
}
\item `Contact` \eqn{{\mu}^{(\rm con)}}{mu_con}: 
\tabular{rrr}{
        \strong{from \\ to} |\tab \strong{adult}               |\tab \strong{child}          \cr
        \strong{adult}      |\tab \emph{reference}             |\tab \code{cont_adult_child} \cr
        \strong{child}      |\tab \code{cont_child_adult} |\tab \code{cont_child_child}
}
}
\if{html}{\out{<hr>}}
}
\section{Augmented Variables}{
\itemize{
\item \code{onset_date} : time of symptom onset of infectious individual;
\item \code{incub_period} : time between infection and symptom onset.
}
\if{html}{\out{<hr>}}
}
\section{Parameters}{
\itemize{
\item \code{beta} : person to person transmission rate;
\item \code{gamma} : contact dilution parameter;
\item \code{sus_vacc} : susceptibility of vaccinated individuals;
\item \code{inf_vacc} : infectivity of vaccinated individuals;
\item \code{cont_adult_child} (c) : contact rate from adult to child;
\item \code{cont_child_adult} (c) : contact rate from child to adult;
\item \code{cont_child_child} (c) : contact rate from adult to child;
\item \code{size_ref} (c) : the size of the reference household;
\item \code{study_period} (c) : study period after inclusion;
\item \code{child_age} (c) : children/adult age threshold.
}
\if{html}{\out{<hr>}}
}
\section{References}{
1. Cauchemez S, Carrat F, Viboud C, Valleron AJ, Boelle PY.
A Bayesian MCMC approach to study transmission of influenza: application to household longitudinal data.
Statist Med. 2004;23: 3469-3487.
\url{https://doi.org/10.1002/sim.1912}
}