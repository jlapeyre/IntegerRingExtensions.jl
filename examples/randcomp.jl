import Random: randstring

function rand_gates(n)
    randstring(['X', 'H', 'S', 'T'], n)
end
