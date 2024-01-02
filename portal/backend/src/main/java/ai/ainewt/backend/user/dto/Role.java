package ai.ainewt.backend.user.dto;

import ai.ainewt.backend.common.EnumRole;
import lombok.Data;
import org.bson.types.ObjectId;

@Data
public class Role {

    private ObjectId id;
    private EnumRole name;

}
