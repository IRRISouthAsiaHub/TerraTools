## -----------------------------------------------------------------------------
library(micsr)


## -----------------------------------------------------------------------------
#| warning: false
trips_2s <- escount(trips + car ~ workschl + size + dist + smsa +
                        fulltime + distnod + realinc + weekend +
                        car | . - car - weekend + adults,
                    data = trips)
names(trips_2s)


## -----------------------------------------------------------------------------
summary(trips_2s)


## -----------------------------------------------------------------------------
#| label: tripsest
#| warning: false
trips_pois <- glm(trips ~ workschl + size + dist + smsa + fulltime +
                      distnod + realinc + weekend + car,
                  data = trips, family = poisson)
trips_ml <- update(trips_2s, method = "ml")


## -----------------------------------------------------------------------------
#| label: drinksest
#| warning: false
#| cache: true
kt_pois <- glm(drinks ~ advice + income + age + educ + race + marital +
                   empstatus + region, data = drinks, family = poisson)
kt_ml <- escount(drinks + advice ~ advice + income + age + educ + race +
                     marital + empstatus + region | income + age + educ +
                     race + marital + empstatus + region + medicare +
                     medicaid + champus + hlthins + regmed + dri + limits +
                     diabete + hearthcond + stroke,
                 data = drinks, method = "ml")
kt_2s <- update(kt_ml, method = "twostep")

