package ai.ainewt.backend.user.dto;

import java.util.Date;
import java.util.HashSet;
import java.util.Set;

import lombok.Data;
import lombok.NoArgsConstructor;
import io.swagger.v3.oas.annotations.media.Schema;
import org.bson.types.ObjectId;

@Data
@NoArgsConstructor
public class User {

    @Schema(description = "bson ObjectId", example = "ObjectId('6580f55ee3315b52c1fd9d2a')", required = false)
    ObjectId id;

    @Schema(description = "User ID", example = "8594125", required = true)
    Long user_id;
    @Schema(description = "User Email Name", example = "some@domain.com", required = true)
    String user_email;
    @Schema(description = "User Name or Nick Name", example = "Jhon Doe", required = true)
    String user_name;
    @Schema(description = "Login Name", example = "JhonDoe", required = false)
    String login_name;
    @Schema(description = "OAuth Provider", example = "Just Email/Google/Apple/Kakao/etc", required = false)
    String oauth_by;
    @Schema(description = "Boolean", example = "True or False", required = false)
    Boolean email_verified;
    @Schema(description = "User Password", example = "SomePassword", required = true)
    String user_password;
    @Schema(description = "First Sign Up Date", example = "ISODate(\"2022-10-05T01:38:08.958Z\")", required = false)
    Date signup_date;
    @Schema(description = "Last Sign In Date", example = "ISODate(\"2022-10-05T01:38:08.958Z\")", required = false)
    Date signin_date;

    private Set<Role> roles = new HashSet<>();

}
