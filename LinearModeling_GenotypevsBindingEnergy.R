library(qvalue)
library(lmtest)
data <- read.csv("EVERYTHING.txt", sep=" ")
attach(data)

sink(file="ROutput.txt")

lm_25 <- lm(log(Kd) ~ x)
print("Protein First Order")
print(coefficients(lm_25))
print(summary(lm_25)$coefficients)
print(summary(lm_25))

lm_26 <- lm(log(Kd) ~ y)
print("Protein site 26")
print(coefficients(lm_26))
print(summary(lm_26)$coefficients)
print(summary(lm_26))

lm_29 <- lm(log(Kd) ~ z)
print("Protein site 29")
print(coefficients(lm_29))
print(summary(lm_29)$coefficients)
print(summary(lm_29))

lm_first_prot <- lm(log(Kd) ~ x + y + z)
print("Protein First Order")
print(coefficients(lm_first_prot))
print(summary(lm_first_prot)$coefficients)
print(summary(lm_first_prot))

lm_first_protDNA <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2)
print("Protein and DNA First Order")
print(coefficients(lm_first_protDNA))
print(summary(lm_first_protDNA)$coefficients)
print(summary(lm_first_protDNA))

xy <- x*y
xz <- x*z
yz <- y*z

lm_second_25_26 <- lm(log(Kd) ~ x + y + z + xy)
print("Protein Second Order 25_26")
print(coefficients(lm_second_25_26))
print(summary(lm_second_25_26)$coefficients)
print(summary(lm_second_25_26))

print(lrtest(lm_second_25_26, lm_first_prot))

lm_second_25_29 <- lm(log(Kd) ~ x + y + z + xz)
print("Protein Second Order 25_29")
print(coefficients(lm_second_25_29))
print(summary(lm_second_25_29)$coefficients)
print(summary(lm_second_25_29))

print(lrtest(lm_second_25_29, lm_first_prot))

lm_second_26_29 <- lm(log(Kd) ~ x + y + z + yz)
print("Protein Second Order 26_29")
print(coefficients(lm_second_26_29))
print(summary(lm_second_26_29)$coefficients)
print(summary(lm_second_26_29))

print(lrtest(lm_second_26_29, lm_first_prot))

lm_second_prot <- lm(log(Kd) ~ x + y + z + xy + xz + yz)
print("Protein Second Order")
print(coefficients(lm_second_prot))
print(summary(lm_second_prot)$coefficients)
print(summary(lm_second_prot))

print(lrtest(lm_second_prot, lm_first_prot))

xyz <- x*y*z

lm_third_prot <- lm(log(Kd) ~ x + y + z + xy + xz + yz + xyz)
print("Protein Third Order")
print(coefficients(lm_third_prot))
print(summary(lm_third_prot)$coefficients)
print(summary(lm_third_prot))

print(lrtest(lm_third_prot, lm_second_prot))

w1w2 <- w1*w2
w1y2 <- w1*y2
w1k2 <- w1*k2
y1w2 <- y1*w2
y1y2 <- y1*y2
y1k2 <- y1*k2
k1w2 <- k1*w2
k1y2 <- k1*y2
k1k2 <- k1*k2

lm_second_protDNA_intra <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2)
print("Protein and DNA Second Order - intra-molecular")
print(coefficients(lm_second_protDNA_intra))
print(summary(lm_second_protDNA_intra)$coefficients)
print(summary(lm_second_protDNA_intra))

print(lrtest(lm_second_protDNA_intra, lm_second_prot))

xw1 <- x*w1
xy1 <- x*y1
xk1 <- x*k1
xw2 <- x*w2
xy2 <- x*y2
xk2 <- x*k2

lm_second_inter_25 <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + xw1 + xy1 + xk1 + xw2 + xy2 + xk2)
print("Inter 25 Order")
print(coefficients(lm_second_inter_25))
print(summary(lm_second_inter_25)$coefficients)
print(summary(lm_second_inter_25))

print(lrtest(lm_second_inter_25, lm_second_protDNA_intra))

yw1 <- y*w1
yy1 <- y*y1
yk1 <- y*k1
yw2 <- y*w2
yy2 <- y*y2
yk2 <- y*k2

lm_second_inter_26 <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + yw1 + yy1 + yk1 + yw2 + yy2 + yk2)
print("Inter 26 Order")
print(coefficients(lm_second_inter_26))
print(summary(lm_second_inter_26)$coefficients)
print(summary(lm_second_inter_26))

print(lrtest(lm_second_inter_26, lm_second_protDNA_intra))

zw1 <- z*w1
zy1 <- z*y1
zk1 <- z*k1
zw2 <- z*w2
zy2 <- z*y2
zk2 <- z*k2

lm_second_inter_29 <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + zw1 + zy1 + zk1 + zw2 + zy2 + zk2)
print("Inter 29 Order")
print(coefficients(lm_second_inter_29))
print(summary(lm_second_inter_29)$coefficients)
print(summary(lm_second_inter_29))

print(lrtest(lm_second_inter_29, lm_second_protDNA_intra))

lm_second_protDNA_inter <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + xw1 + xy1 + xk1 + xw2 + xy2 + xk2 + yw1 + yy1 + yk1 + yw2 + yy2 + yk2 + zw1 + zy1 + zk1 + zw2 + zy2 + zk2)
print("Protein and DNA Second Order - inter-molecular")
print(coefficients(lm_second_protDNA_inter))
print(summary(lm_second_protDNA_inter)$coefficients)
print(summary(lm_second_protDNA_inter))

print(lrtest(lm_second_protDNA_inter, lm_second_protDNA_intra))

xyw1 <- x*y*w1
xyy1 <- x*y*y1
xyk1 <- x*y*k1
xyw2 <- x*y*w2
xyy2 <- x*y*y2
xyk2 <- x*y*k2

lm_third_25_26 <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + xw1 + xy1 + xk1 + xw2 + xy2 + xk2 + yw1 + yy1 + yk1 + yw2 + yy2 + yk2 + zw1 + zy1 + zk1 + zw2 + zy2 + zk2 + xyw1 + xyy1 + xyk1 + xyw2 + xyy2 + xyk2 )
print("Third Order 25_26 - inter-molecular")
print(coefficients(lm_third_25_26))
print(summary(lm_third_25_26)$coefficients)
print(summary(lm_third_25_26))

print(lrtest(lm_third_25_26, lm_second_protDNA_inter))

xzw1 <- x*z*w1
xzy1 <- x*z*y1
xzk1 <- x*z*k1
xzw2 <- x*z*w2
xzy2 <- x*z*y2
xzk2 <- x*z*k2

lm_third_25_29 <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + xw1 + xy1 + xk1 + xw2 + xy2 + xk2 + yw1 + yy1 + yk1 + yw2 + yy2 + yk2 + zw1 + zy1 + zk1 + zw2 + zy2 + zk2 + xzw1 + xzy1 + xzk1 + xzw2 + xzy2 + xzk2 )
print("Third Order 25_29 - inter-molecular")
print(coefficients(lm_third_25_29))
print(summary(lm_third_25_29)$coefficients)
print(summary(lm_third_25_29))

print(lrtest(lm_third_25_29, lm_second_protDNA_inter))

yzw1 <- y*z*w1
yzy1 <- y*z*y1
yzk1 <- y*z*k1
yzw2 <- y*z*w2
yzy2 <- y*z*y2
yzk2 <- y*z*k2

lm_third_26_29 <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + xw1 + xy1 + xk1 + xw2 + xy2 + xk2 + yw1 + yy1 + yk1 + yw2 + yy2 + yk2 + zw1 + zy1 + zk1 + zw2 + zy2 + zk2 + yzw1 + yzy1 + yzk1 + yzw2 + yzy2 + yzk2 )
print("Third Order 26_29 - inter-molecular")
print(coefficients(lm_third_26_29))
print(summary(lm_third_26_29)$coefficients)
print(summary(lm_third_26_29))

print(lrtest(lm_third_26_29, lm_second_protDNA_inter))

xw1w2 <- x*w1*w2
xw1y2 <- x*w1*y2
xw1k2 <- x*w1*k2
xy1w2 <- x*y1*w2
xy1y2 <- x*y1*y2
xy1k2 <- x*y1*k2
xk1w2 <- x*k1*w2
xk1y2 <- x*k1*y2
xk1k2 <- x*k1*k2

lm_third_25_re <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + xw1 + xy1 + xk1 + xw2 + xy2 + xk2 + yw1 + yy1 + yk1 + yw2 + yy2 + yk2 + zw1 + zy1 + zk1 + zw2 + zy2 + zk2 + xw1w2 + xw1y2 + xw1k2 + xy1w2 + xy1y2 + xy1k2 + xk1w2 + xk1y2 + xk1k2 )
print("Third Order 25_re - inter-molecular")
print(coefficients(lm_third_25_re))
print(summary(lm_third_25_re)$coefficients)
print(summary(lm_third_25_re))

print(lrtest(lm_third_25_re, lm_second_protDNA_inter))

yw1w2 <- y*w1*w2
yw1y2 <- y*w1*y2
yw1k2 <- y*w1*k2
yy1w2 <- y*y1*w2
yy1y2 <- y*y1*y2
yy1k2 <- y*y1*k2
yk1w2 <- y*k1*w2
yk1y2 <- y*k1*y2
yk1k2 <- y*k1*k2

lm_third_26_re <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + xw1 + xy1 + xk1 + xw2 + xy2 + xk2 + yw1 + yy1 + yk1 + yw2 + yy2 + yk2 + zw1 + zy1 + zk1 + zw2 + zy2 + zk2 + yw1w2 + yw1y2 + yw1k2 + yy1w2 + yy1y2 + yy1k2 + yk1w2 + yk1y2 + yk1k2 )
print("Third Order 26_re - inter-molecular")
print(coefficients(lm_third_26_re))
print(summary(lm_third_26_re)$coefficients)
print(summary(lm_third_26_re))

print(lrtest(lm_third_26_re, lm_second_protDNA_inter))

zw1w2 <- z*w1*w2
zw1y2 <- z*w1*y2
zw1k2 <- z*w1*k2
zy1w2 <- z*y1*w2
zy1y2 <- z*y1*y2
zy1k2 <- z*y1*k2
zk1w2 <- z*k1*w2
zk1y2 <- z*k1*y2
zk1k2 <- z*k1*k2

lm_third_29_re <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + xw1 + xy1 + xk1 + xw2 + xy2 + xk2 + yw1 + yy1 + yk1 + yw2 + yy2 + yk2 + zw1 + zy1 + zk1 + zw2 + zy2 + zk2 + zw1w2 + zw1y2 + zw1k2 + zy1w2 + zy1y2 + zy1k2 + zk1w2 + zk1y2 + zk1k2 )
print("Third Order 29_re - inter-molecular")
print(coefficients(lm_third_29_re))
print(summary(lm_third_29_re)$coefficients)
print(summary(lm_third_29_re))

print(lrtest(lm_third_29_re, lm_second_protDNA_inter))

lm_third_protDNA_inter <- lm(log(Kd) ~ x + y + z + w1 + y1 + k1 + w2 + y2 + k2 + xy + xz + yz + w1w2 + w1y2 + w1k2 + y1w2 + y1y2 + y1k2 + k1w2 + k1y2 + k1k2 + xw1 + xy1 + xk1 + xw2 + xy2 + xk2 + yw1 + yy1 + yk1 + yw2 + yy2 + yk2 + zw1 + zy1 + zk1 + zw2 + zy2 + zk2 + xyw1 + xyy1 + xyk1 + xyw2 + xyy2 + xyk2 + xzw1 + xzy1 + xzk1 + xzw2 + xzy2 + xzk2 + yzw1 + yzy1 + yzk1 + yzw2 + yzy2 + yzk2 + xw1w2 + xw1y2 + xw1k2 + xy1w2 + xy1y2 + xy1k2 + xk1w2 + xk1y2 + xk1k2 + yw1w2 + yw1y2 + yw1k2 + yy1w2 + yy1y2 + yy1k2 + yk1w2 + yk1y2 + yk1k2 + zw1w2 + zw1y2 + zw1k2 + zy1w2 + zy1y2 + zy1k2 + zk1w2 + zk1y2 + zk1k2)
print("Protein and DNA Third Order - inter-molecular")
print(coefficients(lm_third_protDNA_inter))
print(summary(lm_third_protDNA_inter)$coefficients)
print(summary(lm_third_protDNA_inter))

print(lrtest(lm_third_protDNA_inter, lm_second_protDNA_inter))
