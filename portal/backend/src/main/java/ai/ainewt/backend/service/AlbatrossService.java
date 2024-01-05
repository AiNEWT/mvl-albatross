package ai.ainewt.backend.service;

import org.bson.Document;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

@Service
public class AlbatrossService {
    @Autowired
    AlbatrossDao albatrossDao;

    public Document getContainer(String userName) {
        return albatrossDao.getContainer(userName);
    }

    public Integer getMaxPortNum() {
        return albatrossDao.getMaxPortNum();
    }

    public boolean updateContainerStatus(Document containerStatusDoc) {
        return albatrossDao.updateContainerStatus(containerStatusDoc);
    }

}
