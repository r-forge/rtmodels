"p.error" <-
function (v, a, z = a/2, s2 = 0.1^2) 
{
    p = (exp(-2 * v * a/s2) - exp(-2 * v * z/s2))/(exp(-2 * v * 
        a/s2) - 1)
    p[v == 0] = z/(a - z)
    p
}
