#!/usr/bin/ruby

# ./prof2vtk.rb hoge.prof -> hoge.vtu

def prof2vtk(fname1, fname2)
	prof = File.open(fname1, "r")
	num = prof.gets.to_i
	lines = prof.readlines

	type = []
	pos = []
	vel = []
	press = []
	pressave = []
	i = 0
	for l in lines do
		dat = l.split
		type[i]= dat[1]
		pos[i*3]	= dat[2].to_f
		pos[i*3+1]	= dat[3].to_f
		pos[i*3+2]	= dat[4].to_f
		vel[i*3]	= dat[5].to_f
		vel[i*3+1]	= dat[6].to_f
		vel[i*3+2]	= dat[7].to_f
		press[i]	= dat[8].to_f
		pressave[i]	= dat[9].to_f
		i = i + 1
	end

	puts "generating #{fname2}..."

	vtk = File.open(fname2, "w")
	vtk.puts "<?xml version='1.0' encoding='UTF-8'?>"
	vtk.puts "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>"
	vtk.puts "<UnstructuredGrid>"
	vtk.puts "<Piece NumberOfCells='#{num}' NumberOfPoints='#{num}'>"
	vtk.puts "<Points>"
		vtk.puts "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>"
		for i in 0..num-1 do
			vtk.puts "#{pos[i*3]} #{pos[i*3+1]} #{pos[i*3+2]}"
		end
		vtk.puts "</DataArray>"
	vtk.puts "</Points>"
	vtk.puts "<PointData>"
		vtk.puts "<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>"
			for i in 0..num-1 do
				vtk.puts type[i]
			end
		vtk.puts "</DataArray>"
		vtk.puts "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>"
			for i in 0..num-1 do
				v = Math.sqrt(vel[i*3]*vel[i*3]+ vel[i*3+1]*vel[i*3+1] + vel[i*3+2]*vel[i*3+2])
				vtk.puts "#{v}"
			end
		vtk.puts "</DataArray>"
		vtk.puts "<DataArray NumberOfComponents='1' type='Float32' Name='Press' format='ascii'>"
			for i in 0..num-1 do
				vtk.puts "#{press[i]}"
			end
		vtk.puts "</DataArray>"
		vtk.puts "<DataArray NumberOfComponents='1' type='Float32' Name='PressAverage' format='ascii'>"
			for i in 0..num-1 do
				vtk.puts "#{pressave[i]}"
			end
		vtk.puts "</DataArray>"
	vtk.puts "</PointData>"
	vtk.puts "<Cells>"
		vtk.puts "<DataArray type='Int32' Name='connectivity' format='ascii'>"
			for i in 0..num-1
				vtk.puts "#{i}"
			end
		vtk.puts "</DataArray>"
		vtk.puts "<DataArray type='Int32' Name='offsets' format='ascii'>"
			for i in 0..num-1
				vtk.puts "#{i+1}"
			end
		vtk.puts "</DataArray>"
		vtk.puts "<DataArray type='UInt8' Name='types' format='ascii'>"
			for i in 0..num-1
				vtk.puts "1"
			end
		vtk.puts "</DataArray>"
	vtk.puts "</Cells>"
	vtk.puts "</Piece>"
	vtk.puts "</UnstructuredGrid>"
	vtk.puts "</VTKFile>"
end

if ARGV.size() == 1
	prof2vtk ARGV[0], ARGV[0].gsub(/.prof/, ".vtu")
end
