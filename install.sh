
#ln -s environments/the_counterstream.env.xml

echo "################################"
echo "INSTALLING ENVIRONMENTS & ROBOTS"
echo "################################"
for f in environments/*; do
        d="../../src/data/"`basename $f`
        rm -rf $d
        ln -s $f $d
        echo "ln -s $f $d"
done;
echo "################################"
for f in robots/*; do
        d="../../src/robots/"`basename $f`
        rm -rf $d
        ln -s $f $d
        echo "ln -s $f $d"
done;
echo "################################"
echo "Done"



