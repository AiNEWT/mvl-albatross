package ai.ainewt.backend.user.dto.biz;

import ai.ainewt.backend.auth.JwtResponse;
import ai.ainewt.backend.auth.JwtUtils;
import ai.ainewt.backend.auth.UserDetailsImpl;
import ai.ainewt.backend.common.UserManager;
import ai.ainewt.backend.user.dto.User;
import io.swagger.v3.oas.annotations.Operation;
import io.swagger.v3.oas.annotations.Parameter;
import io.swagger.v3.oas.annotations.Parameters;
import io.swagger.v3.oas.annotations.media.Schema;
import org.bson.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import io.swagger.v3.oas.annotations.tags.Tag;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.security.authentication.AuthenticationManager;
import org.springframework.security.authentication.UsernamePasswordAuthenticationToken;
import org.springframework.security.core.Authentication;
import org.springframework.security.core.context.SecurityContextHolder;
import org.springframework.web.bind.annotation.*;

import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

@Tag(name = "User Controller", description = "Controller for User")
@RestController
public class UserController {
    private static final Logger logger = LoggerFactory.getLogger(UserController.class);

    @Autowired
    AuthenticationManager authenticationManager;
    @Autowired
    JwtUtils jwtUtils;
    @Autowired
    UserService userService;

    @Operation(summary = "User Register", description = "Manage the registration of first-time users")
    @Parameters(value = {
            @Parameter(name = "bodyDoc",required = true,schema = @Schema(implementation = User.class))
        }
    )
    @RequestMapping(value = "/api/v1/user/register",
            consumes = "application/json",
            produces = "application/json",
            method = RequestMethod.POST)
    @ResponseBody
    public ResponseEntity<Document> register(@RequestHeader HttpHeaders headers,
                                             @RequestParam HashMap<String,String> paramsMap,
                                             @RequestBody Document bodyDoc ) throws InterruptedException, ExecutionException {
        ResponseEntity<Document> responseEntity = null;

        boolean checkRequiredfields = false;
        StringBuilder checkMessage = new StringBuilder();

        // Check Inputs
        try {

            if ( bodyDoc.containsKey("username") || bodyDoc.containsKey("email") || bodyDoc.containsKey("password")  ) {
                Document tmp = new Document();
                tmp.append("user_name",bodyDoc.get("username"));
                tmp.append("user_email",bodyDoc.get("email"));
                tmp.append("user_password",bodyDoc.get("password"));
                bodyDoc = tmp;
            }

            logger.debug("bodyDoc : {}", bodyDoc.toString());

            String[] keysToCheck = {"user_email", "user_name", "user_password"};

            for ( String key : keysToCheck ) {
                if (bodyDoc.containsKey(key) && bodyDoc.get(key) != null) {
                    checkRequiredfields = true;
                } else {
                    checkRequiredfields = false;
                    if ( checkMessage.isEmpty() ) {
                        checkMessage.append(key);
                    } else {
                        checkMessage.append(",");
                        checkMessage.append(key);
                    }
                }
            }

        } catch (Exception e) {
            responseEntity = new ResponseEntity<>(
                    new Document("results", e.getMessage()), HttpStatus.INTERNAL_SERVER_ERROR
            );
            return responseEntity;
        }

        if( checkRequiredfields ) {

            boolean checkUserName = userService.checkExistsUserName(bodyDoc.get("user_name", String.class));
            boolean checkUserEmail = userService.checkExistsUserEmail(bodyDoc.get("user_email", String.class));

            if (checkUserName || checkUserEmail) {
                checkMessage.append("User ");
                if (checkUserName)
                    checkMessage.append("Name ");
                if (checkUserEmail)
                    checkMessage.append("Email ");
                checkMessage.append("Already Exists");
                responseEntity = new ResponseEntity<>(new Document("results", checkMessage),HttpStatus.INTERNAL_SERVER_ERROR);
                return responseEntity;
            }

            checkUserEmail = UserManager.isValidEmail(bodyDoc.get("user_email", String.class));
            if(!checkUserEmail) {
                responseEntity = new ResponseEntity<>(
                        new Document("results", String.format("%s is not valid",bodyDoc.get("user_email", String.class))),
                        HttpStatus.INTERNAL_SERVER_ERROR);
                return responseEntity;
            }

            int passswordLength = bodyDoc.get("user_password", String.class).length();

            if(passswordLength < 8) {
                responseEntity = new ResponseEntity<>(
                        new Document("results", "Password More then 8 Characters"),
                        HttpStatus.INTERNAL_SERVER_ERROR);
                return responseEntity;
            }

            Document results = userService.registerUser(bodyDoc);
            responseEntity = new ResponseEntity<>(results,HttpStatus.OK);
        } else {
            logger.debug("checkMessage : {}", checkMessage);
            responseEntity = new ResponseEntity<>(new Document("results", checkMessage),HttpStatus.INTERNAL_SERVER_ERROR);
        }

        return responseEntity;
    }

    @Operation(summary = "User Authenticate", description = "User Authenticate by User Name and Password")
    @Parameters(value = {
            @Parameter(name = "bodyDoc",required = true,schema = @Schema(implementation = User.class))
    }
    )
    @RequestMapping(value = "/api/v1/user/authenticate",
            consumes = "application/json",
            produces = "application/json",
            method = RequestMethod.POST)
    @ResponseBody
    public ResponseEntity<?> authenticate(@RequestHeader HttpHeaders headers,
                                             @RequestParam HashMap<String,String> paramsMap,
                                             @RequestBody Document bodyDoc ) throws InterruptedException, ExecutionException {
        ResponseEntity<?> responseEntity = null;

        boolean checkRequiredfields = false;
        StringBuilder checkMessage = new StringBuilder();

        // Check Inputs
        try {
            logger.debug("bodyDoc : {}", bodyDoc.toString());
            if ( bodyDoc.containsKey("username") || bodyDoc.containsKey("password") ) {
                Document tmp = new Document();
                tmp.append("user_name",bodyDoc.get("username"));
                tmp.append("user_password",bodyDoc.get("password"));
                bodyDoc = tmp;
            }

            String[] keysToCheck = { "user_name", "user_password"};
            for ( String key : keysToCheck ) {
                if (bodyDoc.containsKey(key) && bodyDoc.get(key) != null) {
                    checkRequiredfields = true;
                } else {
                    checkRequiredfields = false;
                    if ( checkMessage.isEmpty() ) {
                        checkMessage.append(key);
                    } else {
                        checkMessage.append(",");
                        checkMessage.append(key);
                    }
                }
            }

            if (!checkRequiredfields) {
                responseEntity = new ResponseEntity<>(new Document("results", checkMessage), HttpStatus.INTERNAL_SERVER_ERROR);
                return responseEntity;
            }

        } catch (Exception e) {
            responseEntity = new ResponseEntity<>(
                    new Document("results", e.getMessage()), HttpStatus.INTERNAL_SERVER_ERROR
            );
            return responseEntity;
        }

        Authentication authentication = authenticationManager.authenticate(
                new UsernamePasswordAuthenticationToken(bodyDoc.get("user_name", String.class), bodyDoc.get("user_password", String.class))
        );

        // TODO : Create a shareable strategy to apply in the future
        SecurityContextHolder.getContext().setAuthentication(authentication);

        String jwt = jwtUtils.generateJwtToken(authentication);
        String jwtRefresh = jwtUtils.generateRefreshJwtToken();

        UserDetailsImpl userDetails = (UserDetailsImpl) authentication.getPrincipal();
        List<String> roles = userDetails.getAuthorities().stream()
                .map(item -> item.getAuthority())
                .collect(Collectors.toList());

        //public JwtResponse(String accessToken, String id, String username, String email, List<String> roles) {
        responseEntity = new ResponseEntity<>(
                new JwtResponse(jwt,jwtRefresh,userDetails.getId(), userDetails.getUsername(), userDetails.getEmail(), roles), HttpStatus.OK
        );

        return responseEntity;
    }

    // TODO : Refresh Token Function
    /*
    @Operation(summary = "Refresh Token", description = "AccessToken Refreshed by Refresh Token")
    @Parameters(value = {
            @Parameter(name = "bodyDoc",required = true,schema = @Schema(implementation = JwtResponse.class))
    }
    )
    @RequestMapping(value = "/api/v1/user/refreshtoken",
            consumes = "application/json",
            produces = "application/json",
            method = RequestMethod.POST)
    @ResponseBody
    public ResponseEntity<?> refreshtoken(@RequestHeader HttpHeaders headers,
                                          @RequestParam HashMap<String,String> paramsMap,
                                          @RequestBody Document bodyDoc ) throws InterruptedException, ExecutionException {
        String requestRefreshToken = bodyDoc.get("refreshToken",String.class);
    }
    */


}
