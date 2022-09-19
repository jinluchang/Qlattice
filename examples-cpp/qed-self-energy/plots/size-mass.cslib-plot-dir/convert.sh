for i in *.mp ; do
fn=${i%.mp}
rm mpost-job.* 2>&1
TEX=latex mpost -jobname mpost-job $fn.mp
for i in mpost-job.? ; do
echo "$i"
mv "$i" "$fn"-"${i#mpost-job.}".eps
done
done
