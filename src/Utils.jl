

function remap(
    nold::Int,
    xold::Array{Float64,1},
    fold::Array{Float64,1},
    nnew::Int,
    xnew::Array{Float64,1},
    fnew::Array{Float64,1},
    av::Int,
    elow::Int,
    eup::Int,
)
    ierror = 0

    if xold[nold] >= xold[1]
        sgnold = 1
        ioldstart = 1
        ioldend = nold
    else
        sgnold = -1
        ioldstart = nold
        ioldend = 1
    end
    if xnew[nnew] >= xnew[1]
        sgnnew = 1
        inewstart = 1
        inewend = nnew
    else
        sgnnew = -1
        inewstart = nnew
        inewend = 1
    end

    for inew = inewstart:inewend:sgnnew
        if xnew[inew] < xold[ioldstart]
            if elow == 0
                ierror = -1
                return ierror
            elseif elow == 1
                fnew[inew] = fold[ioldstart]
            elseif elow == 2
                if nold > 1
                    if ((fold[ioldstart+sgnold] > 0.0) && (fold[ioldstart] > 0.0))
                        fnew[inew] = fold[ioldstart] * (
                            fold[ioldstart+sgnold] / fold[ioldstart]
                        )^(
                            (log(xnew[inew]) - log(xold[ioldstart])) / (log(xold[ioldstart+sgnold]) - log(xold[ioldstart]))
                        )
                    else
                        fnew[inew] = 0.0
                    end
                else
                    fnew[inew] = fold[ioldstart]
                end
            elseif elow == 3
                fnew[inew] = 0.0
            elseif elow == 4
                if nold > 1
                    eps = ((xnew[inew] - xold[ioldstart]) /
                           (xold[ioldstart+sgnold] - xold[ioldstart]))
                    if eps > -1.0
                        fnew[inew] = (1 - abs(eps)) * fold[ioldstart]
                    else
                        fnew[inew] = 0.0
                    end
                else
                    fnew[inew] = fold[ioldstart]
                end
            else
                throw(ArgumentError("Unknown elow"))
            end
        elseif xnew[inew] > xold[ioldend]
            if eup == 0
                ierror = -1
                return ierror
            elseif eup == 1
                fnew[inew] = fold[ioldend]
            elseif eup == 2
                if nold > 1
                    if ((fold[ioldend-sgnold] > 0.0) && (fold[ioldend] > 0.0))
                        fnew[inew] = fold[ioldend] * (
                            fold[ioldend-sgnold] / fold[ioldend]
                        )^(
                            (log(xnew[inew]) - log(xold[ioldend])) / (log(xold[ioldend-sgnold]) - log(xold[ioldend]))
                        )
                    else
                        fnew[inew] = 0.0
                    end
                else
                    fnew[inew] = fold[ioldend]
                end
            elseif eup == 3
                fnew[inew] = 0.0
            elseif eup == 4
                if nold > 1
                    eps = ((xnew[inew] - xold[ioldend]) /
                           (xold[ioldend-sgnold] - xold[ioldend]))
                    if eps > -1.0
                        fnew[inew] = (1 - abs(eps)) * fold[ioldend]
                    else
                        fnew[inew] = 0.0
                    end
                else
                    fnew[inew] = fold[ioldend]
                end
            else
                throw(ArgumentError("Unknown eup"))
            end
        else
            if av == 0
                hunt()
                if xnew[inew] == xold[1]
                    fnew[inew] = fold[1]
                elseif xnew[inew] == xold[nold]
                    fnew[inew] = fold[nold]
                else
                    if ((iold <= 0) || (iold >= nold))
                        eps = (xnew[inew] - xold[iold]) / (xold[iold+1] - xold[iold])
                        if ((fold[iold] > 0.0) && (fold[iold+1] > 0.0))
                            fnew[inew] = exp((1.0 - eps) * log(fold[iold]) + eps * log(fold[iold+1]))
                        else
                            fnew[inew] = (1.0 - eps) * fold[iold] + eps * fold[iold+1]
                        end
                    end
                end
            elseif av == 1
                throw(ArgumentError("Integral conserving mapping not implemented"))
            else
                throw(ArgumentError("Unknown av"))
            end
        end
    end
end