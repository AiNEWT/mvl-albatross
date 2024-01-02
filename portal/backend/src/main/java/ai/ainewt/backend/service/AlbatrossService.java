package ai.ainewt.backend.service;

import ai.ainewt.backend.auth.JwtResponse;
import io.swagger.v3.oas.annotations.Operation;
import io.swagger.v3.oas.annotations.Parameter;
import io.swagger.v3.oas.annotations.Parameters;
import io.swagger.v3.oas.annotations.media.Schema;
import io.swagger.v3.oas.annotations.tags.Tag;
import org.bson.Document;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.security.access.prepost.PreAuthorize;
import org.springframework.web.bind.annotation.*;

import java.util.HashMap;
import java.util.concurrent.ExecutionException;

@Tag(name = "Albatross Controller", description = "Controller for Albatross Docker")
@RestController
public class AlbatrossService {

    @Operation(summary = "Basic User Services", description = "Get User Services")
    @Parameters(value = {
            @Parameter(name = "bodyDoc",required = false,schema = @Schema(implementation = Document.class))
    }
    )
    @RequestMapping(value = "/api/v1/service/user",
//            consumes = "application/json",
            produces = "application/json",
            method = RequestMethod.GET)
    @ResponseBody
    public ResponseEntity<?> getUserService(@RequestHeader HttpHeaders headers
//                                            @RequestParam HashMap<String,String> paramsMap,
//                                            @RequestBody Document bodyDoc
    ) throws InterruptedException, ExecutionException {

        ResponseEntity<?> responseEntity = new ResponseEntity<>(
                "http://localhost:3838/", HttpStatus.OK
        );

        return responseEntity;
    }
}
