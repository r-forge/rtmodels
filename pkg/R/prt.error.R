"prt.error" <-
function (t, v, a, z, s = 0.1, k = 1:100) 
{
    defect = sapply(t, function(t) pi * s^2 * (exp(-v * z/s^2)/a^2) * 
        sum(2 * k * sin(k * pi * z/a) * exp(-(v^2/s^2 + pi^2 * 
            k^2 * s^2/a^2) * t/2)/(v^2/s^2 + pi^2 * k^2 * s^2/a^2)))
    1 - defect/p.error(v = v, a = a, z = z, s2 = s^2)
}
