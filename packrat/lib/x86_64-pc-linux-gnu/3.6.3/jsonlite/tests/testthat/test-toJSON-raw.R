context("toJSON raw")

test_that("Encoding raw vector", {
  x <- list(myraw = charToRaw("bla"))
  x$mydf <- data.frame(foo=1:3)
  x$mydf$bar <- as.character.hexmode(charToRaw("bla"))

  y <- fromJSON(toJSON(x))
  expect_that(x$mydf$bar, is_identical_to(y$mydf$bar))
  expect_that(y$myraw, is_identical_to("Ymxh"))

  # Serialize raw as int
  y <- fromJSON(toJSON(x, raw = 'int'))
  expect_equal(y$myraw, as.integer(x$myraw))

  # Serialize raw as hex
  y <- fromJSON(toJSON(x, raw = 'hex'))
  expect_equal(y$myraw, as.character.hexmode(x$myraw))

  # Serialize raw as JavaScript
  x <- list(myraw = charToRaw("bla"))
  expect_equal(toJSON(x, raw = 'js'), '{"myraw":(new Uint8Array([98,108,97]))}')

})
