package ai.ainewt.backend.service;

import ai.ainewt.backend.user.dto.biz.UserController;
import com.mongodb.client.MongoClient;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import lombok.extern.slf4j.Slf4j;
import org.bson.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Slf4j
@Component
public class AlbatrossDao {

    private static final Logger logger = LoggerFactory.getLogger(AlbatrossDao.class);
    @Autowired
    MongoClient mongoClient;

    public Document getContainer(String userName) {
        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<Document> containerCollection = database.getCollection("containers");

        Document filter = new Document().append("user_name",userName);

        return containerCollection.find(filter).first();
    }

}
