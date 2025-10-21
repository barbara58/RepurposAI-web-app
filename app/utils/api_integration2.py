import requests
import json
import pandas as pd

base_url = "https://api.platform.opentargets.org/api/v4/graphql"

def fetch_primary_drug(strg_chemblid):
  query_string="""
  query MechanismsOfActionSectionQuery($chemblId: String!) {
    drug(chemblId: $chemblId) {
      id
      mechanismsOfAction {
        rows {
          mechanismOfAction
          targetName
          targets {
            id
            approvedSymbol
          }
          references {
            source
            urls
          }
        }
        uniqueActionTypes
        uniqueTargetTypes
      }
      parentMolecule {
        id
        name
      }
      childMolecules {
        id
        name
      }
    }
  }
  """
  variables = {"chemblId": strg_chemblid}

  r = requests.post(base_url, json={"query": query_string, "variables": variables})
  #print(r.status_code)

# Transform API response from JSON into Python dictionary and print in console
  api_response = json.loads(r.text)
 # print(api_response)
  drginfo=[]
  chemblid_res=api_response['data']['drug']['id']
  drginfo.append(chemblid_res)

  for i in api['data']['drug']['mechanismsOfAction']['rows']:
    drginfo.append(i['targets'][0]['id'])
    drginfo.append(i['mechanismOfAction'])
    drginfo.append(i['targetName'])
  
  drginfo_df=pd.DataFrame(drginfo, columns=['chemblid', 'ensemblid', 'mechanism', 'target_name' ])
  return drginfo_df

def fetch_drug_target_associations(strg_ensemblid):
    query_string="""
    query TargetAssociationsQuery(
    $id: String!
    $index: Int!
    $size: Int!
    $sortBy: String!
    $enableIndirect: Boolean!
    $datasources: [DatasourceSettingsInput!]
    $rowsFilter: [String!]
    $facetFilters: [String!]
    $entitySearch: String!
    ) {
    target(ensemblId: $id) {
        id
        approvedSymbol
        associatedDiseases(
        page: { index: $index, size: $size }
        orderByScore: $sortBy
        enableIndirect: $enableIndirect
        datasources: $datasources
        Bs: $rowsFilter
        facetFilters: $facetFilters
        BFilter: $entitySearch
        ) {
        count
        rows {
            disease {
            id
            name
            }
            score
            datasourceScores {
            componentId: id
            score
            }
        }
        }
    }
    }
    """
    variables = {
    "id": strg_ensemblid,
    "index":0,
    "size":300,

    "sortBy": "score",
    "enableIndirect": False,
    "datasources": [
    {
        "id": "gwas_credible_sets",
        "weight": 1,
        "propagate": True,
        "required": False
    },
        {
      "id": "gene_burden",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "eva",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "genomics_england",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "gene2phenotype",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "uniprot_literature",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "uniprot_variants",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "orphanet",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "clingen",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "cancer_gene_census",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "intogen",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "eva_somatic",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "cancer_biomarkers",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "chembl",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "crispr_screen",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "crispr",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "slapenrich",
      "weight": 0.5,
      "propagate": True,
      "required": False
    },
    {
      "id": "progeny",
      "weight": 0.5,
      "propagate": True,
      "required": False
    },
    {
      "id": "reactome",
      "weight": 1,
      "propagate": True,
      "required": False
    },
    {
      "id": "sysbio",
      "weight": 0.5,
      "propagate": True,
      "required": False
    },
    {
      "id": "europepmc",
      "weight": 0.2,
      "propagate": True,
      "required": False
    },
    {
      "id": "expression_atlas",
      "weight": 0.2,
      "propagate": True,
      "required": False
    },
    {
      "id": "impc",
      "weight": 0.2,
      "propagate": True,
      "required": False
    },
    {
      "id": "ot_crispr_validation",
      "weight": 0.5,
      "propagate": True,
      "required": False
    },
    {
      "id": "ot_crispr",
      "weight": 0.5,
      "propagate": True,
      "required": False
    },
    {
      "id": "encore",
      "weight": 0.5,
      "propagate": True,
      "required": False
    }

    ],
    "entity": "target",
    "entitySearch": ""
    }

    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    api_response = json.loads(r.text)
    return api_response

#comment this later
#import sys
def main():
    df_original_molec=fetch_primary_drug(strg_chemblid)
    ensembl_target=df_original_molec['ensemblid']
    alldiseases=fetch_drug_target_associations(ensembl_target)
    disease_data=alldiseases['data']['target']['associatedDiseases']['rows']
    firstlist=list(map(lambda x: [x['disease']['id'],x['disease']['name'],x['score'] ], disease_data))
    #Dataframe with global scores
    df1=pd.DataFrame(firstlist, columns=['id', 'name', 'global_score'])
    allcfds=[]
    for elem in disease_data:
      idd=elem['disease']['id']
      allscores=elem['datasourceScores']
      cdf=pd.DataFrame(allscores).T
      cdf.columns=cdf.iloc[0]
      cdf=cdf.iloc[1:]
      cdf.reset_index(drop=True,inplace=True)
      cdf['id']=idd
      allcfds.append(cdf)
    allotherscores=pd.concat(allcfds, axis=0)
    return df1
main()
