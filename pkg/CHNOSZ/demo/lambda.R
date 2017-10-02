# plot effects of lambda transition in quartz
# after Berman 1988 Figs. 1 and 2
layout(matrix(c(1, 4:2, 1, 7:5), nrow=4), heights=c(0.7, 3, 3, 3))
# plot title first
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "Effects of lambda transition in quartz, after Berman (1988) Figs. 1 and 2", cex=1.8)
opar <- par(mar=c(4, 4.5, 1, 0.5), cex=0.8)

T <- convert(seq(0, 1400, 1), "K")

labplot <- function(x) label.plot(x, xfrac=0.9, yfrac=0.1, paren=TRUE)
# this sets the units used for making the axis labels
E.units("J")
Cplab <- axis.label("Cp")
Vlab <- axis.label("V")
Tlab <- axis.label("T")

# calculate properties at 1 kbar with and without transition
Qz_1bar <- berman("quartz", T=T, units="J")
Qz_1bar_notrans <- berman("quartz", T=T, calc.transition=FALSE, units="J")
# Fig. 1a: volume
plot(T, Qz_1bar$V, type="l", xlab=Tlab, ylab=Vlab, ylim=c(22.5, 24))
legend("topleft", legend="1 bar", bty="n")
labplot("a")
# TODO: why don't we get the curvature his plot for V shows?
# Should it be in the v4 parameter (but it's zero)??

# Fig. 1b: heat capacity
plot(T, Qz_1bar$Cp, type="l", xlab=Tlab, ylab=Cplab)
lines(T, Qz_1bar_notrans$Cp, lty=3)
legend("topleft", legend="1 bar", bty="n")
labplot("b")

# calculate properties at 10 kbar with and without transition
Qz_10bar <- berman("quartz", T=T, P=10000, units="J")
Qz_10bar_notrans <- berman("quartz", T=T, P=10000, calc.transition=FALSE, units="J")
# Fig. 1c: heat capacity
plot(T, Qz_10bar$Cp, type="l", xlab=Tlab, ylab=Cplab)
lines(T, Qz_10bar_notrans$Cp, lty=3)
legend("topleft", legend="10 kb", bty="n")
labplot("c")

# like Ber88 Fig. 2
Tlambda <- 848 # Kelvin
dTdP <- 0.0237
Pkb <- seq(1, 50, 1)
P <- 1000 * Pkb
T <- Tlambda + dTdP * (P - 1)
Qz_withtrans <- berman("quartz", T=T, P=P, units="J")
Qz_notrans <- berman("quartz", T=T, P=P, calc.transition=FALSE, units="J")
Qz_lambda <- Qz_withtrans - Qz_notrans
Plab <- expression(list(italic(P), "kb"))
plot(Pkb, Qz_lambda$G, type="l", ylim=c(-300, -50), ylab=axis.label("DlG"), xlab=Plab)
labplot("d")
plot(Pkb, Qz_lambda$H, type="l", ylim=c(1200, 1800), ylab=axis.label("DlH"), xlab=Plab)
labplot("e")
plot(Pkb, Qz_lambda$S, type="l", ylim=c(0, 3), ylab=axis.label("DlS"), xlab=Plab)
labplot("f")

par(opar)
