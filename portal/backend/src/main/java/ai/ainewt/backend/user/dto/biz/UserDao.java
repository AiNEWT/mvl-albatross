package ai.ainewt.backend.user.dto.biz;

import ai.ainewt.backend.common.EnumRole;
import ai.ainewt.backend.user.dto.Role;
import ai.ainewt.backend.user.dto.User;
import com.mongodb.client.MongoClient;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import com.mongodb.client.model.FindOneAndUpdateOptions;
import com.mongodb.client.model.ReturnDocument;
import com.mongodb.client.result.InsertOneResult;
import lombok.extern.slf4j.Slf4j;
import org.bson.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import java.util.Date;
import java.util.HashSet;
import java.util.Set;

@Slf4j
@Component
public class UserDao {

    @Autowired
    RoleDao roleDao;
    @Autowired
    MongoClient mongoClient;

    // TODO : Apply Transactional
    public Document registerUser(Document userDoc) {
        Document resultsDoc = new Document();
        try {
            MongoDatabase database = mongoClient.getDatabase("portal");
            MongoCollection<Document> userCollection = database.getCollection("user");
            Set<Role> roles = new HashSet<>();

            if (userDoc.containsKey("roles")) {
                String roleName = userDoc.get("role",String.class);
                roles = roleDao.findByName(roleName);
            } else {
                roles = roleDao.findByName("ROLE_USER");
            }

            userDoc.append("signup_date", new Date());
            userDoc.append("roles", roles);
            InsertOneResult insertOneResult = userCollection.insertOne(userDoc);

            resultsDoc.append("results", "OK");
            resultsDoc.append("successMessage", insertOneResult.toString());
        } catch (Exception e) {
            resultsDoc.append("results", "error");
            resultsDoc.append("errorMessage", e.getMessage());
        }
        return resultsDoc;
    }
    public User changePassword(User user) {
        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<User> collection = database.getCollection("user", User.class);
        Document filter = new Document()
                .append("user_id",user.getUser_id())
                .append("user_name",user.getUser_name());

        FindOneAndUpdateOptions findOneAndUpdateOptions =
                new FindOneAndUpdateOptions().upsert(true).returnDocument(ReturnDocument.AFTER);

        return collection.findOneAndUpdate(
                filter,
                new Document("set", new Document("user_password",user.getUser_password())),
                findOneAndUpdateOptions
        );
    }


    public Long getUserId() {
        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<Document> collection = database.getCollection("env");

        FindOneAndUpdateOptions findOneAndUpdateOptions =
                new FindOneAndUpdateOptions().upsert(true).returnDocument(ReturnDocument.AFTER);

        Long addCount = 1L;
        Document userCountDoc = collection.findOneAndUpdate(
          new Document("prop_key","user_count"),
          new Document("$inc", new Document("value",addCount)),
                findOneAndUpdateOptions
        );
        return userCountDoc.get("value",Long.class);
    }

    // TODO : Apply Transactional
    public boolean checkExistsUserId(long checkUserId) {
        boolean exists = false;
        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<Document> collection = database.getCollection("user");
        Document filter = new Document("user_id",checkUserId);
        Document results = collection.find(filter).first();
        if( results != null && results.get("user_id").equals(checkUserId) )  {
            exists = true;
        }
        return exists;
    }

    public User findByUsername(String userName) {
        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<User> collection = database.getCollection("user", User.class);
        Document filter = new Document("user_name",userName);
        return collection.find(filter).first();
    }

    public boolean checkExistsUserName(String userName) {
        boolean exists = false;
        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<Document> collection = database.getCollection("user");
        Document filter = new Document("user_name",userName);
        Document results = collection.find(filter).first();
        if( results != null && results.get("user_name").equals(userName) )  {
            exists = true;
        }
        return exists;
    }

    public boolean checkExistsUserEmail(String userEmail) {
        boolean exists = false;
        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<Document> collection = database.getCollection("user");
        Document filter = new Document("user_email",userEmail);
        Document results = collection.find(filter).first();
        if( results != null && results.get("user_name").equals(userEmail) )  {
            exists = true;
        }
        return exists;
    }

}
