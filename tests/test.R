library(tmda)
print('hello')
# pxkcd(1/sqrt(2*pi)) == 1

# dxkcd(1/sqrt(2*pi)) == 0

# pxkcd(qxkcd(0.1)) - 0.1 < 1e-10

# pxkcd(qxkcd(0.1, swap.end.points = T), swap.end.points = T) - 0.1 < 1e-10

# pxkcd(qxkcd(-0.1, log.p = T), log.p = T) + 0.1 < 1e-10

# pxkcd(qxkcd(-0.1, swap.end.points = T, log.p = T), swap.end.points = T, log.p = T) + 0.1 < 1e-10

# qxkcd(pxkcd(0.1)) - 0.1 < 1e-10

# qxkcd(pxkcd(0.1, swap.end.points = T), swap.end.points = T) - 0.1 < 1e-10

# qxkcd(pxkcd(0.1, log.p = T), log.p = T) - 0.1 < 1e-10

# qxkcd(pxkcd(0.1, swap.end.points = T, log.p = T), swap.end.points = T, log.p = T) - 0.1 < 1e-10
