# stan.fit.ind.20  <- stan.fit
# stan.fit.ind.50  <- stan.fit
# stan.fit.ind.100 <- stan.fit

stan.fit <- stan.fit.ind.100

stan.fit.summary <- summary(stan.fit)[[1]]

est.100 <- stan.fit.summary[grep("bd_delta_eps", dimnames(stan.fit.summary)[[1]]), ][, c(4, 6, 8)] %>%
  as.data.frame() %>% mutate(ind = rownames(.))
names(est.100)[1:3] <- c("lwr", "mid", "upr")
est.100 <- est.100 %>% mutate(num_samples = "100")

stan.fit <- stan.fit.ind.50

stan.fit.summary <- summary(stan.fit)[[1]]

est.50 <- stan.fit.summary[grep("bd_delta_eps", dimnames(stan.fit.summary)[[1]]), ][, c(4, 6, 8)] %>%
  as.data.frame() %>% mutate(ind = rownames(.))
names(est.50)[1:3] <- c("lwr", "mid", "upr")
est.50 <- est.50 %>% mutate(num_samples = "50")

stan.fit <- stan.fit.ind.20

stan.fit.summary <- summary(stan.fit)[[1]]

est.20 <- stan.fit.summary[grep("bd_delta_eps", dimnames(stan.fit.summary)[[1]]), ][, c(4, 6, 8)] %>%
  as.data.frame() %>% mutate(ind = rownames(.))
names(est.20)[1:3] <- c("lwr", "mid", "upr")
est.20 <- est.20 %>% mutate(num_samples = "20")

est.check <- rbind(est.20, est.50, est.100)

est.check %>% group_by(num_samples) %>% mutate(wid = upr - lwr) %>%
  summarize(
    wid_str = sd(wid)
  )

ggplot(est.check, aes(mid, ind, colour = num_samples)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = lwr, xmax = upr))

## Also to confirm, the width of the CI should decrease with the more bd measured

est.20 <- est.20 %>% mutate(
  num_measured = rowSums(bd.measured)
)

est.20 <- est.20 %>% mutate(wid = upr - lwr)

est.20 %>% {
  ggplot(., aes(num_measured, wid)) + geom_point()
}



stan.fit <- stan.fit.1

stan.fit.summary <- summary(stan.fit)[[1]]

est.1 <- stan.fit.summary[
  # grep("bd_delta_eps", dimnames(stan.fit.summary)[[1]])
  1:7, ][, c(4, 6, 8)] %>%
  as.data.frame() %>% mutate(ind = rownames(.))
names(est.1)[1:3] <- c("lwr", "mid", "upr")
est.1 <- est.1 %>% mutate(model = "1")

stan.fit <- stan.fit.2

stan.fit.summary <- summary(stan.fit)[[1]]

est.2 <- stan.fit.summary[
  # grep("bd_delta_eps", dimnames(stan.fit.summary)[[1]])
  1:7, ][, c(4, 6, 8)] %>%
  as.data.frame() %>% mutate(ind = rownames(.))
names(est.2)[1:3] <- c("lwr", "mid", "upr")
est.2 <- est.2 %>% mutate(model = "2")

stan.fit <- stan.fit.3

stan.fit.summary <- summary(stan.fit)[[1]]

est.3 <- stan.fit.summary[
  # grep("bd_delta_eps", dimnames(stan.fit.summary)[[1]])
  1:7, ][, c(4, 6, 8)] %>%
  as.data.frame() %>% mutate(ind = rownames(.))
names(est.3)[1:3] <- c("lwr", "mid", "upr")
est.3 <- est.3 %>% mutate(model = "3")

stan.fit <- stan.fit.4

stan.fit.summary <- summary(stan.fit)[[1]]

est.4 <- stan.fit.summary[
  # grep("bd_delta_eps", dimnames(stan.fit.summary)[[1]])
  1:7, ][, c(4, 6, 8)] %>%
  as.data.frame() %>% mutate(ind = rownames(.))
names(est.4)[1:3] <- c("lwr", "mid", "upr")
est.4 <- est.4 %>% mutate(model = "4")

stan.fit <- stan.fit.6

stan.fit.summary <- summary(stan.fit)[[1]]

est.6 <- stan.fit.summary[
  # grep("bd_delta_eps", dimnames(stan.fit.summary)[[1]])
  1:7, ][, c(4, 6, 8)] %>%
  as.data.frame() %>% mutate(ind = rownames(.))
names(est.6)[1:3] <- c("lwr", "mid", "upr")
est.6 <- est.6 %>% mutate(model = "6")

est.check <- rbind(est.1, est.2, est.3, est.4, est.6)

ggplot(est.check, aes(mid, ind, colour = model)) + 
  geom_point(position = position_dodge(1.2)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), position = position_dodge(1.2), width = 0.2)

