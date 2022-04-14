source etc/bashrc

./Allwmake > logMake.out
sed -i 's/user/'$USER'/g' logMake.out
#######################################################################################i
### check if test cases differents ####
DIFF=$(diff logMake.out logMakeGold.out) 
echo $DIFF
if [ "$DIFF" != "" ] 
then
    echo "Installation error" 

else

    echo "Installation successful" 
    rm logMake.out
    ### delete test cases ##############################################################
fi
