package ai.ainewt.backend.user.dto.biz;

import ai.ainewt.backend.common.EnumRole;
import ai.ainewt.backend.user.dto.Role;
import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoClient;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import lombok.extern.slf4j.Slf4j;
import org.bson.Document;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

@Slf4j
@Component
public class RoleDao {

    @Autowired
    MongoClient mongoClient;
    Set<Role> findByName(String name) {
        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<Role> collection = database.getCollection("roles", Role.class);
        Document filter = new Document("name",name);
        FindIterable<Role> findRoles = collection.find(filter);
        Set<Role> roles = new HashSet<>();
        log.debug("Collection Roles : {}",findRoles.toString());
        for(Role role:findRoles) {
            roles.add(role);
        }
        return roles;
    }
}
