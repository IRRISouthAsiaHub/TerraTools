## -----------------------------------------------------------------------------
#| message: false
library(micsr)
cigmales <- cigmales |>
    transform(age2 = age ^ 2, educ2 = educ ^ 2,
           age3 = age ^ 3, educ3 = educ ^ 3,
           educage = educ * age)                                 
pois_cig <- glm(cigarettes ~ habit + price + restaurant + income + age +
                    age2 + educ + educ2 + famsize + race, data = cigmales,
                family = quasipoisson)


## -----------------------------------------------------------------------------
iv_cig <- expreg(cigarettes ~ habit + price + restaurant + income +
                     age + age2 + educ + educ2 + famsize + race |
                     . - habit + age3 + educ3 + educage + lagprice +
                     reslgth, data = cigmales, method = "iv")
gmm_cig <- update(iv_cig, method = "gmm")


## -----------------------------------------------------------------------------
#| collapse: true
sargan(iv_cig) |> gaze()
sargan(gmm_cig) |> gaze()


## -----------------------------------------------------------------------------
ml_bwt <- glm(birthwt ~ cigarettes + parity + race + sex, data = birthwt,
              family = quasipoisson)
iv_bwt <- expreg(birthwt ~ cigarettes + parity + race + sex |
                     . - cigarettes + edmother + edfather + faminc +
                     cigtax, data = birthwt, method = "iv")
gmm_bwt <- update(iv_bwt, method = "gmm")


## -----------------------------------------------------------------------------
#| collapse: true
sargan(gmm_bwt) |> gaze()

