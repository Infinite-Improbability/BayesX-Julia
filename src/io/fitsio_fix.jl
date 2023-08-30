"""
We need to patch FITSIO because it throws an exception due to the status column on the Chandra data.
This is a really rough fix that just skips the column instead of raising an exception so I'm
not going to contribute it upstream.

The MIT License (MIT)
Copyright (c) 2012 - 2015 FITSIO.jl authors

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
"""

# FITSIO.jl doesn't support BitArray columns and so refuses to read any FITS file containing them.
# Since we have FITS files that do contain them, from Chandra, and we don't care about those columns
# we override the method to bypass this problem. It is a hack, but it works.

@eval FITSIO begin
    ## Arrr me hearties, we're doing some type piracy
    function fits_get_col_info(f::FITSIO.FITSFile, colnum::Integer)
        eqtypecode, repeat, width = FITSIO.fits_get_eqcoltype(f, colnum)
        isvariable = eqtypecode < 0
        eqtypecode = abs(eqtypecode)

        (eqtypecode == 1) && @warn ("BitArray ('X') columns not yet supported. Column $colnum will not be readable.")

        T = CFITSIO_COLTYPE_var[eqtypecode]

        if isvariable
            if T !== String
                T = Vector{T}
            end
            rowsize = Int[]
        else
            if T === String
                # for strings, cfitsio only considers it to be a vector column if
                # width != repeat, even if tdim is multi-valued.
                if repeat == width
                    rowsize = Int[]
                else
                    tdim = FITSIO.fits_read_tdim(f, colnum)
                    # if tdim isn't multi-valued, ignore it (we know it *is* a
                    # vector column). If it is multi-valued, prefer it to repeat
                    # width.
                    if length(tdim) == 1
                        rowsize = [div(repeat, width)]
                    else
                        rowsize = tdim[2:end]
                    end
                end
            else
                if repeat == 1
                    rowsize = Int[]
                else
                    rowsize = FITSIO.fits_read_tdim(f, colnum)
                end
            end
        end

        return T, rowsize, isvariable
    end
end