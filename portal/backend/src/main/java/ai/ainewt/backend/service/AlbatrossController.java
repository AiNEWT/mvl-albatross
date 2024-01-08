package ai.ainewt.backend.service;

import io.swagger.v3.oas.annotations.Operation;
import io.swagger.v3.oas.annotations.Parameter;
import io.swagger.v3.oas.annotations.Parameters;
import io.swagger.v3.oas.annotations.media.Schema;
import io.swagger.v3.oas.annotations.tags.Tag;
import org.bson.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.*;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.client.RestTemplate;

import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;

@Tag(name = "Albatross Controller", description = "Controller for Albatross Docker")
@RestController
public class AlbatrossController {

    private static final Logger logger = LoggerFactory.getLogger(AlbatrossController.class);

    @Autowired
    AlbatrossService albatrossService;

    private Document getContainerStatus(String containerId) {
        RestTemplate restTemplate = new RestTemplate();
        String inspectContainerUrl = String.format(
                "http://localhost:2375/containers/%s/json",
                containerId);
        ResponseEntity<String> response
                = restTemplate.getForEntity(inspectContainerUrl, String.class);
        //logger.debug("inspectContainerUrl Document: {}",response);
        Document inspectContainer = Document.parse(response.getBody());
        return  inspectContainer;
    }

    @Operation(summary = "Basic User Services", description = "Get User Services")
    @Parameters(value = {
            @Parameter(name = "bodyDoc",required = false,schema = @Schema(implementation = Document.class))
    }
    )
    @RequestMapping(value = "/api/v1/service/user",// consumes = "application/json",
            produces = "application/json",
            method = RequestMethod.GET)
    @ResponseBody
    public ResponseEntity<?> getUserService(@RequestHeader HttpHeaders headers,
                                            @RequestParam HashMap<String,String> paramsMap ) throws InterruptedException, ExecutionException {

        ResponseEntity<?> responseEntity = null;

        if(paramsMap.containsKey("username")) {

            logger.debug("Username : {}",paramsMap.get("username"));

            try {

                Document userContainer = albatrossService.getContainer(paramsMap.get("username"));
                //logger.debug("userContainer Document: {}",userContainer.toString());
                if( userContainer == null) {
                    responseEntity = new ResponseEntity<>(new Document("Status", "None"), HttpStatus.OK);
                    return responseEntity;
                }
                StringBuilder frameUrl = new StringBuilder();
                frameUrl.append("http://localhost:");
                frameUrl.append(userContainer.get("ports").toString());
                //logger.debug("frameUrl : {}",frameUrl);

                Document inspectContainer = getContainerStatus(userContainer.get("container_id").toString());

                Document containerState = inspectContainer.get("State",Document.class);
                logger.debug("containerState Document: {}",containerState);

                if( containerState.get("Running").equals(true) &&
                        containerState.get("Paused").equals(false) &&
                        containerState.get("Restarting").equals(false) ) {
                    containerState.append("frame_url", frameUrl);
                    responseEntity = new ResponseEntity<>(containerState, HttpStatus.OK);
                } else {
                    containerState.append("frame_url", "");
                    responseEntity = new ResponseEntity<>(containerState, HttpStatus.OK);
                }

            } catch (Exception e) {
                logger.error("getUserService error : {}", e.getMessage());
                responseEntity = new ResponseEntity<>(new Document("Status", "None"), HttpStatus.OK);
            }

        } else {

            responseEntity = new ResponseEntity<>(new Document("Status", "None"), HttpStatus.OK);

        }

        return responseEntity;
    }

    @Operation(summary = "Basic User Services", description = "Create User Services")
    @Parameters(value = {
            @Parameter(name = "bodyDoc",required = false,schema = @Schema(implementation = Document.class))
    }
    )
    @RequestMapping(value = "/api/v1/service/create",
            consumes = "application/json",
            produces = "application/json",
            method = RequestMethod.POST)
    @ResponseBody
    public ResponseEntity<?> createUserService(@RequestHeader HttpHeaders headers,
                                               //@RequestParam HashMap<String,String> paramsMap,
                                               @RequestBody Document bodyDoc ) throws InterruptedException, ExecutionException {

        ResponseEntity<?> responseEntity = null;
        //logger.debug("createContainerDoc bodyDoc : {}",bodyDoc.toString());

        String userName = bodyDoc.getString("username");
        RestTemplate restTemplate = new RestTemplate();

        // Check Existing Cluster
        Document userContainer = albatrossService.getContainer(bodyDoc.get("username").toString());
        if (userContainer != null && !userContainer.isEmpty()) {
            Document inspectContainer = getContainerStatus(userContainer.get("container_id").toString());
            Document containerState = inspectContainer.get("State",Document.class);

            if (containerState.get("Running").equals(false) ) {
                HttpEntity<String> request = new HttpEntity<String>(headers);
                String restartContainerUrl = String.format("http://localhost:2375/containers/%s/start",userContainer.get("container_id").toString());
                HttpEntity<String> restartRequest = new HttpEntity<String>(headers);
                ResponseEntity<Document> responseContainerReStart = restTemplate.postForEntity(restartContainerUrl,restartRequest,Document.class);
                responseEntity = new ResponseEntity<>(userContainer, HttpStatus.OK);
            } else {
                responseEntity = new ResponseEntity<>(userContainer, HttpStatus.OK);
                return responseEntity;
            }
        }

        Integer maxPortNum = albatrossService.getMaxPortNum();

        String createContainerUrl = "http://localhost:2375/containers/create";

        Document createContainerDoc = new Document()
                .append("Image","mvl-albatross:test")
                .append("ExposedPorts",new Document(String.format("%d/tcp",maxPortNum),new Document()))
                .append("HostConfig",
                        new Document("PortBindings",
                                new Document(String.format("%d/tcp",maxPortNum),
                                        Arrays.asList(new Document("HostPort","3838"))
                                )
                            )
                        )
                .append("RestartPolicy", new Document("Name","always"))
                .append("AutoRemove",true)
                ;

        logger.debug("createContainerDoc : {}",createContainerDoc.toString());

        HttpHeaders postHeaders = new HttpHeaders();
        postHeaders.setContentType(MediaType.APPLICATION_JSON);

        try {

            HttpEntity<String> request = new HttpEntity<String>(createContainerDoc.toJson(),headers);

            //logger.debug("createContainerDoc request: {}",request);

            ResponseEntity<String> response = restTemplate.postForEntity(createContainerUrl,request,String.class);
            Document inspectContainer = Document.parse(response.getBody());

            //logger.debug("createContainerDoc inspectContainer: {}",inspectContainer);

            String containerId = inspectContainer.get("Id",String.class);
            String containerIdShort = containerId.substring(0,12);

            Document containerDoc = new Document();
            containerDoc.append("user_name",userName);
            containerDoc.append("container_id",containerIdShort);
            containerDoc.append("ports",maxPortNum);
            //containerDoc.append("name","");

            String startContainerUrl = String.format("http://localhost:2375/containers/%s/start",containerIdShort);

            request = new HttpEntity<String>(headers);
            ResponseEntity<Document> responseContainerStart = restTemplate.postForEntity(startContainerUrl,request,Document.class);

            //logger.debug("createContainerDoc responseContainerStart: {}",responseContainerStart);
            boolean updateContainerStatus = albatrossService.updateContainerStatus(containerDoc);

            responseEntity = new ResponseEntity<>(containerDoc, HttpStatus.OK);
        } catch (Exception e) {
            logger.error("createContainerDoc error : {}", e.getMessage());
            responseEntity = new ResponseEntity<>(new Document("Status", "None"), HttpStatus.INTERNAL_SERVER_ERROR);
        }

        return responseEntity;
    }
}
