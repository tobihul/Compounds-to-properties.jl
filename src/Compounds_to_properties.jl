module Compounds_to_properties

#The following code can be used to get chemical properties found on PubChem for any list of compounds
#Make sure your list of compounds is in a CSV file and the first cell on A1 is: 'Compounds' or something similar so it recognizes that column as the list
#The list can be in the form of: Chemical names or SMILES but make sure to use the right function for it, it cannot be a combination of the two
#To obtain CAS # for your compounds a separate function is needed and a separate CSV will be created
using Statistics, CSV, DataFrames ,LinearAlgebra, PubChemCrawler
import HTTP
const prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
#This function converts the chemical names to PubChem CID's and is implemeted in the next two functions
function my_get_cid(; name=nothing, smiles=nothing, kwargs...)
    input = "compound/"
    name !== nothing && (input *= "name/$(HTTP.escapeuri(name))/")
    smiles !== nothing && (input *= "smiles/$((smiles))/")
    url = prolog * input * "cids/TXT"
    r = HTTP.request("GET", url; kwargs...)
    cids_string = String(r.body)
    cids = split(cids_string, "\n")
    cids = [cid for cid in cids if !isempty(cid) && !isspace(cid[1])]
    return parse(Int, cids[1])
end
function get_properties_chemical_name(Compounds::Matrix{String})
    cids_f::Vector{Int32} = zeros(length(Compounds))
    for i = 1:length(cids)
        cids_f[i] = my_get_cid(name=Compounds[i])
        @show i
    end
    return cids_f
end
function get_properties_smiles(Compounds::Matrix{String})
    cids_f::Vector{Int32} = zeros(length(Compounds))
    for i = 1:length(cids)
        cids_f[i] = my_get_cid(smiles=Compounds[i])
        @show i
    end
    return cids_f
end
#This function is used to convert CID's to CAS #
function cid_to_cas(cid::Int32)
    pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$cid/xrefs/RN/TXT"
    response = HTTP.get(pubchem_url)
    cas_numbers = split(String(response.body), "\n")
    return cas_numbers[1]
end
function cas_gen(Cids::Vector{Int32})
    CAS = Vector{String}(undef, length(Cids))
    for i = 1:length(Cids)
        CAS[i] = cid_to_cas(Cids[i])
        @show i
    end
    return CAS
end


#Load in your data using this example
#Change the file destination to the location where the CSV with the compounds is located
#The following works for both Strings and names"

Compounds = Matrix{String}(CSV.read("/Users/tobias/Library/CloudStorage/OneDrive-UvA/Research project UvA/Pestmix all compounds.csv", DataFrame))

#Next, if your list is chemical names, use the following function to convert them to CID's

Cids = get_properties_chemical_name(Compounds)

#If your compounds are in the form of SMILES use the following instead

Cids = get_properties_smiles(Compounds)

#If you encounter an error, check in the REPL which iteration it was (the compound in question will be on the row number in the csv corresponding to
#the iteration that failed plus 2. So if the error was after i = 2 then go to row 4 of the csv, change the name or smiles to a different synonym or correct 
#the typo/mistake, save the csv and run all lines up to this point it again

#The next step is to get the properties from PubChem out of the CID's into a DataFrame
#This line needs to be run two times, the first time there will be an I/O error and then it will work

#In this line you can also change which properties you want to be returned. Check PubChem for all the properties that are available and type it into
# the part that says properties. Here is an example: = "property1, property2, property3". In this case I return MonoisotopicMass,XlogP and InChIKey.
properties = CSV.File(get_for_cids(Cids; properties="MonoisotopicMass,XlogP,InChIKey", output="CSV")) |> DataFrame

properties  #To view the DataFrame type display(properties) in the REPL

#Use this line to save the CSV at a desired file location

CSV.write("/Users/tobias/Documents/pestproperties all.csv", properties)

#Use the following to get the CAS # for all compounds in a separate CSV, note that this function is very slow (roughly 0.7 seconds per compound)

CAS = cas_gen(Cids)

CAS_df = DataFrame(CAS = CAS)

CSV.write("/Users/tobias/Documents/All cas of compounds.csv", CAS_df)

export my_get_cid, get_properties_chemical_name, get_properties_smiles, cid_to_cas, cas_gen
end
