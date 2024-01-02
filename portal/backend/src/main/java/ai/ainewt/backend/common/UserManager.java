package ai.ainewt.backend.common;

import ai.ainewt.backend.user.dto.biz.UserService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import java.util.UUID;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@Component
public class UserManager {
    @Autowired
    UserService userService;

    public UserManager() {
    }

    // 이메일 주소의 유효성을 검사하는 정규 표현식
    private static final String EMAIL_REGEX =
            "^[a-zA-Z0-9_+&*-]+(?:\\.[a-zA-Z0-9_+&*-]+)*@(?:[a-zA-Z0-9-]+\\.)+[a-zA-Z]{2,7}$";

    private static final Pattern pattern = Pattern.compile(EMAIL_REGEX);

    /**
     * 이메일 주소가 유효한지 검사하는 함수
     * @param email 검사할 이메일 주소
     * @return 이메일 주소가 유효하면 true, 그렇지 않으면 false 반환
     */
    public static boolean isValidEmail(String email) {
        if (email == null) {
            return false;
        }
        Matcher matcher = pattern.matcher(email);
        return matcher.matches();
    }


}
