########################
## Some project notes ##
########################

## Some model fitting considerations:

 ## 1) Model has gotten really slow when fitting all populations. Some ways of dealing with this:
  ## A) Fitting individual species by themselves since species is just fit as a fixed effect anyway
  ## B) Collapsing the sampling of FL to week

## Some potential model changes:

## 1) Could pursue a better model for Bd load (e.g. a zero inflated regression for example), and may want to use unscaled Bd
 ##   -- There are a few positives to this:
  ##    A) scaled Bd makes it really weird to estimate survival when an individual is uninfected, because
   ##      with the current strategy "uninfected" doesn't really have a defined meaning as Bd is simply a continuous state
   ##      and anything greater than 2 sd from the mean in this continuous load isn't very sensible.
  ##    B) Doing a zero-inflated model would lead to more sensible estimates of what "uninfected" means
 ##   -- However, I am skeptical as to the importance of this change for a few reasons:
  ##    A) The goal is to estimate an individuals yearly average or max load -- most individuals do get infected. Having zero's are mostly
   ##      due to the observation window and do not provide a good estimate of true value within the year
  ##    B) With so much observation error, I would be very skeptical on what a zero means
  ##    C) In sum, given the observation error and that most individuals get infected, I think the continuous load only model is reasonably sensible 

## 2) There is probably a better way to measure population sizes given the issue with estimating population size on days with 0 captures
  ##    A) May want to estimate less frequently than every day, maybe population size using average captures at the level of the random effect
   ##     (intersection of primary period and population)?

## 3) Could give some more effort to exploring exactly when the MeHg model can be fit
 ##      Just because a low proportion of individuals get measured doesn't _immediately_ mean it isn't enough to estimate the effect

## 4) Maybe just discussion points? but the difference in the width of the CI for the main effect of Bd between the 
 ##      MeHg interaction and Bd only model makes me a bit nervous 

