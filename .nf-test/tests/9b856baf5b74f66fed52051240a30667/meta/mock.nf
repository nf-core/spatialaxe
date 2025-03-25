import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies

include { UNZIP  } from '/Users/obrovkina/Documents/Projects/spatialxe/modules/nf-core/unzip/main.nf'


// include test process
include { POINTS2REGIONS_CLUSTER } from '/Users/obrovkina/Documents/Projects/spatialxe/modules/local/points2regions/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    
    {
        def input = []
        
                input[0] = [[], file('https://raw.githubusercontent.com/nf-core/test-datasets/spatialxe/Xenium_Prime_Mouse_Ileum_tiny_outs.zip', checkIfExists: true)]
                
        UNZIP(*input)
    }
    

    // process mapping
    def input = []
    
                input[0] = Channel.of([
                        [id: "test"],
                    ]).combine(UNZIP.out.unzipped_archive.map { it[1] }  + "/transcripts.csv")
                input[1] = 20
                input[2] = 5
                
    //----

    //run process
    POINTS2REGIONS_CLUSTER(*input)

    if (POINTS2REGIONS_CLUSTER.output){

        // consumes all named output channels and stores items in a json file
        for (def name in POINTS2REGIONS_CLUSTER.out.getNames()) {
            serializeChannel(name, POINTS2REGIONS_CLUSTER.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = POINTS2REGIONS_CLUSTER.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
