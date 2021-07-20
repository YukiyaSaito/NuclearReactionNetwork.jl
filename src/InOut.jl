module InOut
using ..ReactionTypes




# begin
# 	open("moller2003.dat") do file
# 	index = 1
# 		for ln in eachline(file)
# 			println([parse(Int, ss) for ss in split(ln)])
			
# 			index += 1
# 		end
# 	end
# end


# using PyCall
# ng = fort.FortranFile("input/nuclear/Reaclib/moller2003.bin","r")
# begin
# 	ng.read_ints()
# 	typeof(ng.read_ints())
# end

end