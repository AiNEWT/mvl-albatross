package ai.ainewt.backend.service;

import ai.ainewt.backend.user.dto.biz.UserController;
import com.mongodb.client.AggregateIterable;
import com.mongodb.client.MongoClient;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import com.mongodb.client.model.FindOneAndUpdateOptions;
import com.mongodb.client.model.ReturnDocument;
import com.mongodb.client.model.UpdateOptions;
import com.mongodb.client.result.UpdateResult;
import lombok.extern.slf4j.Slf4j;
import org.bson.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import java.util.Arrays;

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
    public boolean updateContainerStatus(Document containerStatusDoc) {
        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<Document> containerCollection = database.getCollection("containers");

        Document filter = new Document()
                .append("user_name",containerStatusDoc.get("user_name"));

        Document updateSet = new Document()
                .append("$set",
                        new Document()
                                .append("container_id",containerStatusDoc.get("container_id"))
                                .append("ports",containerStatusDoc.get("ports"))
                        );

        UpdateOptions updateOptions = new UpdateOptions().upsert(true);

        UpdateResult updateResult = containerCollection.updateOne(
                filter,
                updateSet,
                updateOptions
        );
        return updateResult.wasAcknowledged();
    }

    public Integer getMaxPortNum() {

        Integer maxPortNum = 3838;

        MongoDatabase database = mongoClient.getDatabase("portal");
        MongoCollection<Document> containerCollection = database.getCollection("containers");

        AggregateIterable<Document> result =
                containerCollection.aggregate(
                        Arrays.asList(new Document("$match",
                        new Document("ports",
                                new Document("$gte", 3838L))),
                        new Document("$group",
                                new Document("_id", "max_port_num")
                                        .append("ports",
                                                new Document("$max", "$ports"))))
                );

        try {
            if( result.iterator().hasNext() ) {
                Document maxDoc = result.iterator().next();
                String portsString = maxDoc.get("ports").toString();
                maxPortNum = Integer.parseInt(portsString) + 1;
            } else {
                maxPortNum = 3838;
            }
        } catch (Exception e) {
            logger.error("getMaxPortNum error : {}",e.getMessage());
        }

        return maxPortNum;

    }

}
