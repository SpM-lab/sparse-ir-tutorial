# Python
pip3 install --upgrade pip
pip3 install -r requirements.txt

# Julia
julia --project=@. -e "using Pkg; Pkg.instantiate()"
