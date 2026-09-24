library(dfidx)
oopts <- options(dfidx.pos_idx = 1L)

# check whether the same results is obtained with `[`:
# - using a character, a numeric or a logical
# - having the index in the first or another position
# - using one or two arguments

# [ with only one argument (columns)

## idx in first position
mn <- dfidx(munnell[1:8], c(region = "state", president = "year"))
d1 <- mn[c(3, 5)]
expect_equal(unname(idx_name(d1)), 1L)
d2 <- mn[c(1, 3, 5)]
d3 <- mn[-c(1, 2, 4)]
d4 <- mn[- c(2, 4)]
d5 <- mn[c(FALSE, FALSE, TRUE, FALSE, TRUE)]
d6 <- mn[c(TRUE, FALSE, TRUE, FALSE, TRUE)]
d7 <- mn[c("idx", "highway", "utilities")]
d8 <- mn[c("highway", "utilities")]
d9 <- subset(mn, select = c("highway", "utilities"))
d10 <- subset(mn, select = c("highway", "utilities", "idx"))
d <- sapply(list(d2, d3, d4, d5, d6, d7, d8, d9, d10), function(x) identical(x, d1))
expect_true(all(d))

## idx in second position
mn <- dfidx(munnell[1:8], c(region = "state", president = "year"), position = 2L, name = "index")
d1 <- mn[c(3, 5)]
expect_equal(unname(idx_name(d1)), 3L)
d2 <- mn[c(2, 3, 5)]
d3 <- mn[-c(1, 2, 4)]
d4 <- mn[- c(1, 4)]
d5 <- mn[c(FALSE, FALSE, TRUE, FALSE, TRUE)]
d6 <- mn[c(FALSE, TRUE, TRUE, FALSE, TRUE)]
d7 <- mn[c("index", "highway", "utilities")]
d8 <- mn[c("highway", "utilities")]
d9 <- subset(mn, select = c("highway", "utilities"))
d10 <- subset(mn, select = c("highway", "utilities", "index"))
d <- sapply(list(d2, d3, d4, d5, d6, d7, d8, d9, d10), function(x) identical(x, d1))
expect_true(all(d))

# [ with two arguments (lines and columns) slice

## idx in first position
mn <- dfidx(munnell[1:8], c(region = "state", president = "year"), name = "index")
d1 <- mn[1:10, c(3, 5)]
expect_equal(unname(idx_name(d1)), 1L)
d2 <- mn[1:10, c(1, 3, 5)]
d3 <- mn[1:10, -c(1, 2, 4)]
d4 <- mn[1:10, - c(2, 4)]
d5 <- mn[1:10, c(FALSE, FALSE, TRUE, FALSE, TRUE)]
d6 <- mn[1:10, c(TRUE, FALSE, TRUE, FALSE, TRUE)]
d7 <- mn[1:10, c("index", "highway", "utilities")]
d8 <- mn[1:10, c("highway", "utilities")]
d9 <- subset(mn, 1:10, select = c("highway", "utilities"), drop.unused.levels = FALSE)
d10 <- subset(mn, 1:10, select = c("highway", "utilities", "index"), drop.unused.levels = FALSE)
d <- sapply(list(d2, d3, d4, d5, d6, d7, d8, d9, d10), function(x) identical(x, d1))
expect_true(all(d))

## idx in second position
mn <- dfidx(munnell[1:8], c(region = "state", president = "year"), position = 2L, name = "index")
d1 <- mn[1:10, c(3, 5)]
expect_equal(unname(idx_name(d1)), 3L)
d2 <- mn[1:10, c(2, 3, 5)]
d3 <- mn[1:10, - c(1, 2, 4)]
d4 <- mn[1:10, - c(1, 4)]
d5 <- mn[1:10, c(FALSE, FALSE, TRUE, FALSE, TRUE)]
d6 <- mn[1:10, c(FALSE, TRUE, TRUE, FALSE, TRUE)]
d7 <- mn[1:10, c("index", "highway", "utilities")]
d8 <- mn[1:10, c("highway", "utilities")]
d9 <- subset(mn, 1:10, select = c("highway", "utilities"), drop.unused.levels = FALSE)
d10 <- subset(mn, 1:10, select = c("highway", "utilities", "index"), drop.unused.levels = FALSE)
d <- sapply(list(d2, d3, d4, d5, d6, d7, d8, d9, d10), function(x) identical(x, d1))
expect_true(all(d))

# [ with two arguments (lines and columns) filter

## idx in first position
mn <- dfidx(munnell[1:8], c(region = "state", president = "year"))
d1 <- mn[mn$water > 4E03, c(3, 5)]
expect_equal(unname(idx_name(d1)), 1L)
d2 <- mn[mn$water > 4E03, c(1, 3, 5)]
d3 <- mn[mn$water > 4E03, -c(1, 2, 4)]
d4 <- mn[mn$water > 4E03, - c(2, 4)]
d5 <- mn[mn$water > 4E03, c(FALSE, FALSE, TRUE, FALSE, TRUE)]
d6 <- mn[mn$water > 4E03, c(TRUE, FALSE, TRUE, FALSE, TRUE)]
d7 <- mn[mn$water > 4E03, c("idx", "highway", "utilities")]
d8 <- mn[mn$water > 4E03, c("highway", "utilities")]
d9 <- subset(mn, water > 4E03, select = c("highway", "utilities"), drop.unused.levels = FALSE)
d10 <- subset(mn, water > 4E03, select = c("highway", "utilities", "idx"), drop.unused.levels = FALSE)
d <- sapply(list(d2, d3, d4, d5, d6, d7, d8, d9, d10), function(x) identical(x, d1))
expect_true(all(d))

## idx in second position
mn <- dfidx(munnell[1:8], c(region = "state", president = "year"), position = 2L, name = "index")
d1 <- mn[mn$water > 4E03, c(3, 5)]
expect_equal(unname(idx_name(d1)), 3L)
d2 <- mn[mn$water > 4E03, c(2, 3, 5)]
d3 <- mn[mn$water > 4E03, - c(1, 2, 4)]
d4 <- mn[mn$water > 4E03, - c(1, 4)]
d5 <- mn[mn$water > 4E03, c(FALSE, FALSE, TRUE, FALSE, TRUE)]
d6 <- mn[mn$water > 4E03, c(FALSE, TRUE, TRUE, FALSE, TRUE)]
d7 <- mn[mn$water > 4E03, c("index", "highway", "utilities")]
d8 <- mn[mn$water > 4E03, c("highway", "utilities")]
d9 <- subset(mn, mn$water > 4E03, select = c("highway", "utilities"), drop.unused.levels = FALSE)
d10 <- subset(mn, mn$water > 4E03, select = c("highway", "utilities", "index"), drop.unused.levels = FALSE)
d <- sapply(list(d2, d3, d4, d5, d6, d7, d8, d9, d10), function(x) identical(x, d1))
expect_true(all(d))

## extracting one series

mn1 <- mn[, "highway", drop = TRUE]
mn2 <- mn[["highway"]]
mn3 <- mn$highway
d <- sapply(list(mn2, mn3), function(x) identical(x, mn1))
expect_true(all(d))

# transform

mn1 <- transform(mn, lwater = log(water), wat77 = ifelse(year == 1977, water, 0))
mn2 <- mn
mn2$lwater <- log(mn2$water)
mn2$wat77 <- (mn2$index$year == 1977) * mn2$water
expect_identical(mn1, mn2)

# organize

mn1 <- organize(mn, year, water)
mn2 <- mn[order(mn$year, mn$water), ]
expect_identical(mn1, mn2)

mn3 <- organize(mn, year, - water)
mn4 <- mn[order(mn$year, - mn$water), ]
expect_identical(mn3, mn4)

options(oopts)
