package ai.ainewt.backend.user.dto.biz;

import ai.ainewt.backend.user.dto.User;
import org.bson.Document;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.crypto.password.PasswordEncoder;
import org.springframework.stereotype.Service;

@Service
public class UserService {
    @Autowired
    UserDao userDao;

    @Autowired
    PasswordEncoder passwordEncoder; // DI

    public Document registerUser(Document UserDoc) {
        //long userId = userManager.getUserId();
        Long userId = getUserId();
        if (!checkExistsUserId(userId)){
            UserDoc.append("user_id",userId);
            UserDoc.replace("user_password", passwordEncoder.encode(UserDoc.get("user_password").toString()));
            return userDao.registerUser(UserDoc);
        } else {
            Document resultsDoc = new Document();
            resultsDoc.append("results", "error");
            resultsDoc.append("errorMessage", "User Id Already Exists");
            return resultsDoc;
        }
    }

    public User changePassword(User user) {
        return userDao.changePassword(user);
    }

    public Long getUserId() {
        return userDao.getUserId();
    }

    public boolean checkExistsUserId(long checkUserId) {
        return userDao.checkExistsUserId(checkUserId);
    }
    public boolean checkExistsUserName(String userName) {
        return userDao.checkExistsUserName(userName);
    }
    public boolean checkExistsUserEmail(String userEmail) {
        return userDao.checkExistsUserEmail(userEmail);
    }

}
