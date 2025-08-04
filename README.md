# What is this?
prediction_crafting.py is used for generating tangents when there are two or more circles representing predictions of typhoon at different times.

There are two functions contributing to the purpose: generator() and pinner().

Pinner() takes two circles of longitude, latitudes and radii and calculate the end of the tangent connecting the circles.

Generator() takes in an array, which contains arrays of length 3 conataining the lon, lat and radii of the circles and pass it through pinner() pair by pair.

Then the function will generate a smooth path done by a cubic spline.

There is a sample data attached inside and matplotlib is used for visualising the results.
