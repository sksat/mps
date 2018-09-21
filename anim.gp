set term gif animate
set output "out.gif"
set size ratio 1
n = 0
n1= system("find ".ARG1." -type f -name '*.prof' | sort | tail -n1 | sed -e 's/[^0-9]//g'")
dn= @ARG2
size = @ARG3

set xrange[-size*0.1:size*0.5]
set yrange[-size*0.1:size*0.5]

while(n<=n1){
	frame = sprintf("%s/output%010d.prof", ARG1, n)
	plot frame every ::1 u 3:5 w p title frame
	n = n + dn
}
