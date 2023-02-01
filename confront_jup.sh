#! /usr/bin/zsh
NAME=$1
echo "mv ${NAME}.py conf_${NAME}.py"
mv ${NAME}.py conf_${NAME}.py
echo "conv_jup ${NAME}.ipynb"
jupyter nbconvert --to script ${NAME}.ipynb
echo "diff ${NAME}.py conf_${NAME}.py"
diff ${NAME}.py conf_${NAME}.py
echo "mv conf_${NAME}.py ${NAME}.py # everything back as before"
mv conf_${NAME}.py ${NAME}.py
